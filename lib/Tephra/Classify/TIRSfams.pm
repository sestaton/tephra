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
     'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::GT';


=head1 NAME

Tephra::Classify::TIRSams - Classify TIR transposons into superfamilies

=head1 VERSION

Version 0.09.8

=cut

our $VERSION = '0.09.8';
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
sub find_tc1_mariner {
    my $self = shift;
    my ($feature, $header, $index, $log) = @_;
    my $fasta = $self->genome->absolute->resolve;
    my $gff   = $self->gff->absolute->resolve;

    my @lengths;
    my $mar_feats;
    my $is_mariner = 0;
    my $has_pdoms  = 0;
    my $pdoms      = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_tc1-mariner.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_tc1-mariner.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_tc1-mariner_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";;
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand, @tirs);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
	my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	for my $tir_feature (@{$feature->{$rep_region}}) {
	    if ($tir_feature->{type} eq 'protein_match') {
		$has_pdoms = 1;
		push @all_pdoms, $tir_feature->{attributes}{name}[0];
	    }

	    if ($tir_feature->{type} eq 'target_site_duplication') {
		($seq_id, $source, $start, $end, $strand) = @{$tir_feature}{qw(seq_id source start end strand)};
		my $tsd_len = $end - $start + 1;
		if ($tsd_len == 2) {
		    my $seq = $self->subseq($index, $seq_id, undef, $start, $end, undef);
		    $is_mariner = 1 if $seq =~ /ta/i;
		}
	    }
	}
	    
	if ($is_mariner) {
	    push @tirs, $rep_region;
	    say $out join "\t", $feature->{$rep_region}[0]{seq_id}, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

	    for my $tir_feature (@{$feature->{$rep_region}}) {
		if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		    my $id = join "_", 'DTT', $elem_id, $seq_id, $start, $end;
		    say $faout join "\n", ">".$id, $seq;

		    $tir_feature->{attributes}{superfamily} = 'DTT';
		    $len = $end - $start + 1;
		}
		my $gff3_str = gff3_format_feature($tir_feature);
		chomp $gff3_str;
		say $out $gff3_str;
	    }

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
    close $out;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) { 
	delete $feature->{$_} for @tirs; 

	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = defined $stat->mean ? sprintf("%.2f", $stat->mean) : 0;
	my $count = $stat->count;

	my $pad = ' ' x 10;
	$log->info("Results - Total number of Tc1-Mariner elements:$pad           $count");
	$log->info("Results - Minimum length of Tc1-Mariner elements:$pad         $min");
	$log->info("Results - Maximum length of Tc1-Mariner elements:$pad         $max");
	$log->info("Results - Mean length of Tc1-Mariner elements:$pad            $mean");
	$log->info("Results - Number of Tc1-Mariner elements protein matches:$pad $pdoms");

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;	
	return (undef, undef, 0);
    }
}

sub find_hat {
    my $self = shift;
    my ($feature, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;
    
    my @lengths;
    my $hat_feats;
    my $is_hat = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_hAT.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_hAT.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_hAT_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand, @tirs);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
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
	    push @tirs, $rep_region;
	    say $out join "\t", $feature->{$rep_region}[0]{seq_id}, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

	    for my $tir_feature (@{$feature->{$rep_region}}) {
                if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		    my $id = join "_", 'DTA', $elem_id, $seq_id, $start, $end;
		    say $faout join "\n", ">".$id, $seq;

		    $tir_feature->{attributes}{superfamily} = 'DTA';
		    $len = $end - $start + 1;
		}
		my $gff3_str = gff3_format_feature($tir_feature);
		chomp $gff3_str;
		say $out $gff3_str;
	    }
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
    close $out;
    close $faout;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;    

    if (@lengths) { 
	delete $feature->{$_} for @tirs; 

	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = defined $stat->mean ? sprintf("%.2f", $stat->mean) : 0;
	my $count = $stat->count;

	my $pad = ' ' x 10;
	$log->info("Results - Total number of hAT elements:$pad                   $count");
        $log->info("Results - Minimum length of hAT elements:$pad                 $min");
        $log->info("Results - Maximum length of hAT elements:$pad                 $max");
        $log->info("Results - Mean length of hAT elements:$pad                    $mean");
	$log->info("Results - Number of hAT elements protein matches:$pad         $pdoms");

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
    }
}

sub find_mutator {
    my $self = shift;
    my ($feature, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;

    my @lengths;
    my $mut_feats;
    my $is_mutator = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_mutator.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_mutator.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_mutator_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand, @tirs);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
            if ($tir_feature->{type} eq 'protein_match') {
                $has_pdoms = 1;
                push @all_pdoms, $tir_feature->{attributes}{name}[0];
            }

            if ($tir_feature->{type} eq 'target_site_duplication') {
                my ($seq_id, $tsd_start, $tsd_end) = @{$tir_feature}{qw(seq_id start end)};
                my $tsd_len = $tsd_end - $tsd_start + 1;
		if ($tsd_len >= 8 && $tsd_len <= 11) {
		    $is_mutator = 1;
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};    
		}
	    }
	}

        if ($is_mutator) {
	    push @tirs, $rep_region;
            say $out join "\t", $feature->{$rep_region}[0]{seq_id}, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

	    for my $tir_feature (@{$feature->{$rep_region}}) {
		if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		    my $id = join "_", 'DTM', $elem_id, $seq_id, $start, $end;
		    say $faout join "\n", ">".$id, $seq;

		    $tir_feature->{attributes}{superfamily} = 'DTM';
		    $len = $end - $start + 1;
		}
		my $gff3_str = gff3_format_feature($tir_feature);
		chomp $gff3_str;
		say $out $gff3_str;
	    }
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
    close $out;
    close $faout;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) { 
	delete $feature->{$_} for @tirs; 

	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = defined $stat->mean ? sprintf("%.2f", $stat->mean) : 0;
	my $count = $stat->count;

	my $pad = ' ' x 10;
	$log->info("Results - Total number of Mutator elements:$pad               $count");
        $log->info("Results - Minimum length of Mutator elements:$pad             $min");
        $log->info("Results - Maximum length of Mutator elements:$pad             $max");
        $log->info("Results - Mean length of Mutator elements:$pad                $mean");
	$log->info("Results - Number of Mutator elements protein matches:$pad     $pdoms");

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
    }
}

sub find_cacta {
    my $self = shift;
    my ($feature, $header, $index, $log) = @_;
    my $fasta = $self->genome->absolute->resolve;
    my $gff   = $self->gff->absolute->resolve;

    my @lengths;
    my $cac_feats;
    my $is_cacta = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_cacta.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_cacta.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_cacta_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand, $tir_len, @tirs);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
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
                my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		if ($seq =~ /^cact(?:a|g)?|cact(?:a|g)?$/i) {
		    # Lewin, 1997 http://www.plantphysiol.org/content/132/1/52.full
		    # provides this TIR length definition, but it seems to remove all predictions
		    $is_cacta = 1; # if $tir_len >= 10 && $tir_len <= 28;
		}		
            }
        }

        if ($is_cacta) {
	    push @tirs, $rep_region;
            say $out join "\t", $feature->{$rep_region}[0]{seq_id}, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";

	    for my $tir_feature (@{$feature->{$rep_region}}) {
		if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		    my $elem_id = $tir_feature->{attributes}{ID}[0];
		    ($seq_id, $source, $start, $end, $strand) = 
			@{$tir_feature}{qw(seq_id source start end strand)};
		    my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		    my $id = join "_", 'DTC', $elem_id, $seq_id, $start, $end;
		    say $faout join "\n", ">".$id, $seq;

		    $tir_feature->{attributes}{superfamily} = 'DTC';
		    $len = $end - $start + 1;
		}
		my $gff3_str = gff3_format_feature($tir_feature);
		chomp $gff3_str;
		say $out $gff3_str;
	    }
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
    close $out;
    close $faout;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) {
	delete $feature->{$_} for @tirs; 

	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = defined $stat->mean ? sprintf("%.2f", $stat->mean) : 0;
	my $count = $stat->count;
	
	my $pad = ' ' x 10;
	$log->info("Results - Total number of CACTA elements:$pad                 $count");
        $log->info("Results - Minimum length of CACTA elements:$pad               $min");
        $log->info("Results - Maximum length of CACTA elements:$pad               $max");
        $log->info("Results - Mean length of CACTA elements:$pad                  $mean");
	$log->info("Results - Number of CACTA elements protein matches:$pad       $pdoms");

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
    }
}

sub write_unclassified_tirs {
    my $self = shift;
    my ($feature, $header, $index, $log) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;

    my @lengths;
    my $unc_feats;
    my $is_unclass = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_unclassified.gff3');
    my $fas        = File::Spec->catfile($path, $name.'_unclassified.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_unclassified_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
            if ($tir_feature->{type} eq 'protein_match') {
                $has_pdoms = 1;
                push @all_pdoms, $tir_feature->{attributes}{name}[0];
            }

            if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
                my $elem_id = $tir_feature->{attributes}{ID}[0];
                ($seq_id, $source, $start, $end, $strand) = 
		    @{$tir_feature}{qw(seq_id source start end strand)};
                my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
                my $id = join "_", 'DTX', $elem_id, $seq_id, $start, $end;
		$tir_feature->{attributes}{superfamily} = 'DTX';

                $lines .= join "\n", ">".$id, $seq;
                $len = $end - $start + 1;
            }
            my $gff3_str = gff3_format_feature($tir_feature);
            $unc_feats .= $gff3_str;
	}
	
	chomp $unc_feats;
	say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
	say $out $unc_feats;
	say $faout $lines;

	undef $unc_feats;
	undef $lines;

	delete $feature->{$rep_region};
	push @lengths, $len;
	$pdom_org = join ",", @all_pdoms;
	$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	$pdoms++ if $has_pdoms;
	$has_pdoms = 0;
    }
    close $out;
    close $faout;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {            
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;
    
    if (@lengths) { 
	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = defined $stat->mean ? sprintf("%.2f", $stat->mean) : 0;
	my $count = $stat->count;

	my $pad = ' ' x 10;
	$log->info("Results - Total number of unclassified TIR elements:$pad      $count");
        $log->info("Results - Minimum length of unclassified TIR elements:$pad    $min");
        $log->info("Results - Maximum length of unclassified TIR elements:$pad    $max");
        $log->info("Results - Mean length of unclassified TIR elements:$pad       $mean");
	$log->info("Results - Number of unclas. TIR elements protein matches:$pad $pdoms");

	return ($outfile, $fas, $count);
    }
    else {
	unlink $outfile, $fas;
	return (undef, undef, 0);
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
	my $lines = do { 
	    local $/ = undef; 
	    open my $fh_in, '<', $file or die "\n[ERROR]: Could not open file: $file\n";
	    <$fh_in>;
	};
	print $out $lines;
    }
    close $out;

    my $gt  = $self->get_gt_exec;
    my $cmd = "$gt gff3 -sort @{$outfiles->{gffs}} > $outfile";
    #say STDERR "debug: $cmd";
    my @out = capture { system([0..5], $cmd) };
    unlink @{$outfiles->{fastas}}, @{$outfiles->{gffs}};

    return;
}

sub subseq {
    my $self = shift;
    my ($index, $loc, $elem, $start, $end, $out) = @_;

    my $location = "$loc:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);
    croak "\n[ERROR]: Something went wrong. This is a bug, please report it.\n"
        unless $length;

    $seq =~ s/.{60}\K/\n/g;

    return $seq;
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
