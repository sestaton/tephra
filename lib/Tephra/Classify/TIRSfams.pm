package Tephra::Classify::TIRSfams;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::GFF3::LowLevel qw(gff3_format_feature);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use Cwd                 qw(getcwd abs_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Path::Class::File;
use Try::Tiny;
use Carp 'croak';
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Logger',
     'Tephra::Role::File',
     'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::GT',
     'Tephra::Classify::Role::LogResults';


=head1 NAME

Tephra::Classify::TIRSams - Classify TIR transposons into superfamilies

=head1 VERSION

Version 0.12.2

=cut

our $VERSION = '0.12.2';
$VERSION = eval $VERSION;

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has gff => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has outfile => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

#
# methods
#
sub find_mites {
    my $self = shift;
    my ($features, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;

    my ($has_pdoms, $mite_num) = (0, 0);
    my (%pdom_index, %mite_index, @mites, @all_pdoms, @mite_lengths, $mite_feats);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $moutfile = File::Spec->catfile($path, $name.'_mite.gff3');
    my $mfas     = File::Spec->catfile($path, $name.'_mite.fasta');

    open my $mout, '>>', $moutfile or die "\n[ERROR]: Could not open file: $moutfile\n";
    open my $mfaout, '>>', $mfas or die "\n[ERROR]: Could not open file: $mfas\n";
    say $mout $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand);
    for my $chr_id (nsort keys %$features) { 
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	    for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($tir_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		}
		
		if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
		    #my $id = join "_", 'DTI', $elem_id, $seq_id, $start, $end;
		    #$tir_feature->{attributes}{superfamily} = 'DTI';
		    my $id = join "_", 'DTX', $elem_id, $seq_id, $start, $end;
		    $tir_feature->{attributes}{superfamily} = 'DTX';
		    
		    $lines .= join "\n", ">".$id, $seq;
		}
		my $gff3_str = gff3_format_feature($tir_feature);
		$mite_feats .= $gff3_str;
	    }
	    chomp $mite_feats;

	    my $len = $end - $start + 1;
	    if (!$has_pdoms && $len <= 600) {
		# MITEs are typically small; hard-coded at 600 bp max, and lack coding domains
		# --
		# Tourist ref: https://www.ncbi.nlm.nih.gov/pubmed/1332797/
		# --
		# Heartbreaker ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC15555/
		# --
		# Feschotte C, Zhang X, Wessler S. Miniature inverted-repeat transposable elements (MITEs) 
		# and their relationship with established DNA transposons. 2002.
		# --
		# Review: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2167627/
		$mite_num++;
		my $unc_mite_feats = $self->_format_mite_features($mite_feats, $lines, $mite_num);
		say $mout join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
		say $mout $unc_mite_feats->{unc_feats};
		say $mfaout $unc_mite_feats->{new_fas};
		
		push @mite_lengths, $len;
		$mite_index{ $unc_mite_feats->{old_id} } = $unc_mite_feats->{new_id};
		push @mites, join "||", $seq_id, $rep_region;
	    }	    
	    undef $mite_feats;
	    undef $lines;

	    $has_pdoms = 0;
	}
    }
    close $mout;
    close $mfaout;
    #dd \@mites;
    
    if (@mite_lengths) {
	$self->remove_repeat_region_features($features, \@mites);

	my $mite_count = 
	    $self->log_basic_element_stats({ lengths => \@mite_lengths, type => 'MITE', log => $log, pdom_ct => 0 });

	return ({ mite_outfile  => $moutfile,
                  mite_fasta    => $mfas,
		  mite_count    => $mite_count, 
		  mite_index    => \%mite_index, });
    }
    else {
	unlink $moutfile, $mfas;

	return undef;
    }
}

sub find_tc1_mariner {
    my $self = shift;
    my ($features, $header, $index, $log) = @_;
    my $fasta = $self->genome->absolute->resolve;
    my $gff   = $self->gff->absolute->resolve;

    my ($is_mariner, $has_pdoms, $pdoms) = (0, 0, 0);
    my (%pdom_index, $pdom_org, $mar_feats, @all_pdoms, @lengths);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_tc1-mariner.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_tc1-mariner.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_tc1-mariner_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    say $out $header;

    my ($lines, $seq_id, $source, $start, $end, $strand, @tirs);
    for my $chr_id (nsort keys %$features) {
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	    for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($tir_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    push @all_pdoms, $tir_feature->{attributes}{name}[0];
		}
		
		if ($tir_feature->{type} eq 'target_site_duplication') {
		    ($seq_id, $source, $start, $end, $strand) = @{$tir_feature}{qw(seq_id source start end strand)};
		    my $tsd_len = $end - $start + 1;
		    if ($tsd_len == 2) {
			my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
			$is_mariner = 1 if $seq =~ /ta/i;
		    }
		}
	    }
	    
	    if ($is_mariner) {
		push @tirs, join "||", $seq_id, $rep_region;
		say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

		for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
			my $elem_id = $tir_feature->{attributes}{ID}[0];
			($seq_id, $source, $start, $end, $strand) = 
			    @{$tir_feature}{qw(seq_id source start end strand)};
			my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
			my $id = join "_", 'DTT', $elem_id, $seq_id, $start, $end;
			say $faout join "\n", ">".$id, $seq;
			
			$tir_feature->{attributes}{superfamily} = 'DTT';
		    }
		    my $gff3_str = gff3_format_feature($tir_feature);
		    chomp $gff3_str;
		    say $out $gff3_str;
		}
		
		my $len = $end - $start + 1;
		push @lengths, $len;

		$pdom_org = join ",", @all_pdoms;
		$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
		$pdoms++ if $has_pdoms;
	    }	    
	    undef $pdom_org;
	    
	    @all_pdoms  = ();
	    $is_mariner = 0;
	    $has_pdoms  = 0;
	}
    }
    close $out;

    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
	if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) {
	$self->remove_repeat_region_features($features, \@tirs);
	my $count = $self->log_basic_element_stats({ lengths => \@lengths, type => 'Tc1-Mariner', log => $log, pdom_ct => $pdoms });

	return ($outfile, $fas, $count); 
    }
    else {
	unlink $outfile, $fas;	
	return (undef, undef, 0);
    }
}

sub find_hat {
    my $self = shift;
    my ($features, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;
    
    my ($is_hat, $has_pdoms, $pdoms) = (0, 0, 0);
    my (%pdom_index, @lengths, @all_pdoms, $hat_feats, $pdom_org);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_hAT.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_hAT.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_hAT_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    say $out $header;

    my ($lines, $seq_id, $source, $start, $end, $strand, @tirs);
    for my $chr_id (nsort keys %$features) {
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	    for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($tir_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    push @all_pdoms, $tir_feature->{attributes}{name}[0];
		}
		
		if ($tir_feature->{type} eq 'target_site_duplication') {
		    ($seq_id, $start, $end, $source, $strand) = 
			@{$tir_feature}{qw(seq_id start end source strand)};
		    my $tsd_len = $end - $start + 1;
		    $is_hat = 1 if $tsd_len == 8;
		}
	    }
	    
	    if ($is_hat) {
		push @tirs, join "||", $seq_id, $rep_region;
		say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

		for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
			my $elem_id = $tir_feature->{attributes}{ID}[0];
			($seq_id, $source, $start, $end, $strand) = 
			    @{$tir_feature}{qw(seq_id source start end strand)};
			my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
			my $id = join "_", 'DTA', $elem_id, $seq_id, $start, $end;
			say $faout join "\n", ">".$id, $seq;
			
			$tir_feature->{attributes}{superfamily} = 'DTA';
		    }
		    my $gff3_str = gff3_format_feature($tir_feature);
		    chomp $gff3_str;
		    say $out $gff3_str;
		}

		my $len = $end - $start + 1;
		push @lengths, $len;

		$pdom_org = join ",", @all_pdoms;
		$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
		$pdoms++ if $has_pdoms;
	    }
	    undef $pdom_org;
	    
	    @all_pdoms = ();
	    $is_hat    = 0;
	    $has_pdoms = 0;
	}
    }
    close $out;
    close $faout;

    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
        if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) { 
	$self->remove_repeat_region_features($features, \@tirs);
	my $count = $self->log_basic_element_stats({ lengths => \@lengths, type => 'hAT', log => $log, pdom_ct => $pdoms });

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
    }
}

sub find_mutator {
    my $self = shift;
    my ($features, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;

    my ($is_mutator, $has_pdoms, $pdoms) = (0, 0, 0);
    my (%pdom_index, @lengths, @all_pdoms, $pdom_org, $mut_feats);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_mutator.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_mutator.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_mutator_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    say $out $header;

    my ($lines, $seq_id, $source, $start, $end, $strand, @tirs);
    for my $chr_id (nsort keys %$features) { 
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	    for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($tir_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    push @all_pdoms, $tir_feature->{attributes}{name}[0];
		}
		
		if ($tir_feature->{type} eq 'target_site_duplication') {
		    my ($tsd_start, $tsd_end) = @{$tir_feature}{qw(start end)};
		    my $tsd_len = $tsd_end - $tsd_start + 1;
		    if ($tsd_len >= 8 && $tsd_len <= 11) {
			$is_mutator = 1;
			($seq_id, $source, $start, $end, $strand) = 
			    @{$tir_feature}{qw(seq_id source start end strand)};    
		    }
		}
	    }
	    
	    if ($is_mutator) {
		push @tirs, join "||", $seq_id, $rep_region;
		say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

		for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
			my $elem_id = $tir_feature->{attributes}{ID}[0];
			($seq_id, $source, $start, $end, $strand) = 
			    @{$tir_feature}{qw(seq_id source start end strand)};
			my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
			my $id = join "_", 'DTM', $elem_id, $seq_id, $start, $end;
			say $faout join "\n", ">".$id, $seq;
			
			$tir_feature->{attributes}{superfamily} = 'DTM';
		    }
		    my $gff3_str = gff3_format_feature($tir_feature);
		    chomp $gff3_str;
		    say $out $gff3_str;
		}

		my $len = $end - $start + 1;
		push @lengths, $len;
		
		$pdom_org = join ",", @all_pdoms;
		$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
		$pdoms++ if $has_pdoms;
	    }
	    undef $pdom_org;
	    
	    @all_pdoms  = ();
	    $is_mutator = 0;
	    $has_pdoms  = 0;
	}
    }
    close $out;
    close $faout;

    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
        if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) { 
	$self->remove_repeat_region_features($features, \@tirs);
	my $count = $self->log_basic_element_stats({ lengths => \@lengths, type => 'Mutator', log => $log, pdom_ct => $pdoms });

	return ($outfile, $fas, $count); 
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
    }
}

sub find_cacta {
    my $self = shift;
    my ($features, $header, $index, $log) = @_;
    my $fasta = $self->genome->absolute->resolve;
    my $gff   = $self->gff->absolute->resolve;

    my ($is_cacta, $has_pdoms, $pdoms) = (0, 0, 0);
    my (%pdom_index, @all_pdoms, @lengths, $cac_feats, $pdom_org);
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_cacta.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_cacta.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_cacta_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand, $tir_len, @tirs);
    for my $chr_id (nsort keys %$features) { 
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	    for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($tir_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    push @all_pdoms, $tir_feature->{attributes}{name}[0];
		}
		
		if ($tir_feature->{type} eq 'terminal_inverted_repeat') {
		    my ($tir_start, $tir_end) = @{$tir_feature}{qw(start end)};                            
		    $tir_len = $tir_end - $tir_start + 1; 
		}
		
		if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
		    
		    if ($seq =~ /^cact(?:a|g)?|cact(?:a|g)?$/i) {
			# Lewin, 1997 http://www.plantphysiol.org/content/132/1/52.full
			# provides this TIR length definition, but it seems to remove all predictions
			$is_cacta = 1; # if $tir_len >= 10 && $tir_len <= 28;
		    }		
		}
	    }

	    if ($is_cacta) {
		push @tirs, join "||", $seq_id, $rep_region;
		say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
		
		for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
			my $elem_id = $tir_feature->{attributes}{ID}[0];
			($seq_id, $source, $start, $end, $strand) = 
			    @{$tir_feature}{qw(seq_id source start end strand)};
			my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
			my $id = join "_", 'DTC', $elem_id, $seq_id, $start, $end;
			say $faout join "\n", ">".$id, $seq;
			
			$tir_feature->{attributes}{superfamily} = 'DTC';
		    }
		    my $gff3_str = gff3_format_feature($tir_feature);
		    chomp $gff3_str;
		    say $out $gff3_str;
		}

		my $len = $end - $start + 1;
		push @lengths, $len;
		
		$pdom_org = join ",", @all_pdoms;
		$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
		$pdoms++ if $has_pdoms;
	    }
	    undef $pdom_org;

	    @all_pdoms = ();
	    $is_cacta  = 0;
	    $has_pdoms = 0;
	}
    }
    close $out;
    close $faout;

    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
        if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) {
	$self->remove_repeat_region_features($features, \@tirs);
	my $count = $self->log_basic_element_stats({ lengths => \@lengths, type => 'CACTA', log => $log, pdom_ct => $pdoms });

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
    }
}

sub write_unclassified_tirs {
    my $self = shift;
    my ($features, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;

    my ($is_unclass, $has_pdoms, $pdoms) = (0, 0, 0);
    my (%pdom_index, @all_pdoms, @unc_lengths, $pdom_org, $unc_feats);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_unclassified.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_unclassified.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_unclassified_domain_org.tsv');

    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand);
    for my $chr_id (nsort keys %$features) { 
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	    for my $tir_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($tir_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    push @all_pdoms, $tir_feature->{attributes}{name}[0];
		}
		
		if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
		    my $id = join "_", 'DTX', $elem_id, $seq_id, $start, $end;
		    $tir_feature->{attributes}{superfamily} = 'DTX';
		    
		    $lines .= join "\n", ">".$id, $seq;
		}
		my $gff3_str = gff3_format_feature($tir_feature);
		$unc_feats .= $gff3_str;
	    }	    
	    chomp $unc_feats;
	    say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

	    say $faout $lines;
	    say $out $unc_feats;
		
	    my $len = $end - $start + 1;
	    push @unc_lengths, $len;
	    
	    undef $unc_feats;
	    undef $lines;

	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	    $has_pdoms = 0;
	}
    }
    close $out;
    close $faout;
    
    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
        if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;
    
    if (@unc_lengths) { 
	my $unc_count = 
	    $self->log_basic_element_stats({ lengths => \@unc_lengths, type => 'unclassified TIR', log => $log, pdom_ct => $pdoms });

	return ({ unc_outfile  => $outfile, 
		  unc_fasta    => $fas, 
		  unc_count    => $unc_count }); 
    }
    else {
	unlink $outfile, $fas;

	return undef;
    }
}

sub write_combined_output {
    my $self = shift;
    my ($outfiles) = @_;
    my $outfile = $self->outfile;

    my ($name, $path, $suffix) = fileparse($outfile, qr/\.[^.]*/);
    my $fasout = File::Spec->catfile($path, $name.'.fasta');
    open my $out, '>', $fasout or die "\n[ERROR]: Could not open file: $fasout\n";

    for my $file (@{$outfiles->{fastas}}) {
	$self->collate($file, $out);
    }
    close $out;

    my $gt  = $self->get_gt_exec;
    my $cmd = "$gt gff3 -sort @{$outfiles->{gffs}} > $outfile";
    #say STDERR "debug: $cmd";
    my @out = capture { system([0..5], $cmd) };
    unlink @{$outfiles->{fastas}}, @{$outfiles->{gffs}};

    return;
}

sub _format_mite_features {
    my $self = shift;
    my ($unc_feats, $lines, $mite_num) = @_;

    my ($new_feats, $old_id);
    my $is_mite = 0;
    my $mite_type = 'MITE';
    my $new_id = $mite_type.$mite_num;

    for my $feat (split /^/, $unc_feats) {
        chomp $feat;
        my @feats = split /\t/, $feat;
	if ($feats[2] eq 'terminal_inverted_repeat_element') {
            $feats[2] = $mite_type;
            $is_mite = 1;
        }
        if ($feats[8] =~ /(terminal_inverted_repeat_element\d+)/) {
	    $old_id = $1;
	    $feats[8] =~ s/$old_id/$new_id/g;  
        }
        $new_feats .= join "\t", @feats, "\n";
    }
    chomp $new_feats;

    my $new_lines;
    for my $line (split /^/, $lines) {
	if ($line =~ />/) { 
	    $line =~ s/$old_id/$new_id/g;
	    $new_lines .= $line."\n";
	}
	else {
	    $new_lines .= $line."\n";
	}
    }
    chomp $new_lines;

    return ({ unc_feats => $new_feats, is_mite => $is_mite, old_id => $old_id, new_id => $new_id, new_fas => $new_lines });
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Classify::TIRSFams


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
