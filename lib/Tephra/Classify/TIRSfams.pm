package Tephra::Classify::TIRSfams;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use IPC::System::Simple qw(capture);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use Path::Class::File;
use Try::Tiny;
use Cwd;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

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

#
# methods
#
sub find_tc1_mariner {
    my $self = shift;
    my ($feature, $header) = @_;
    my $fasta = $self->genome;
    my $gff   = $self->gff;
    
    my @lengths;
    my $mar_feats;
    my $is_mariner = 0;
    my $has_pdoms  = 0;
    my $pdoms      = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_tc1-mariner.gff3");
    open my $out, '>>', $outfile;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    $has_pdoms = 1 if $feats[2] =~ /protein_match/;
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		if ($tsd_len == 2) {
		    #my $tsd = substr $hash->{$loc}, $feats[3], $tsd_len;
		    my $tmp = $tir.".fasta";
		    my $cmd = "samtools faidx $fasta $region:$feats[3]-$feats[4] > $tmp";
		    $self->run_cmd($cmd);
		    my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
		    while (my $seqobj = $seqio->next_seq) {
			my $seq = $seqobj->seq;
			#$seq =~ s/.{60}\K/\n/g;
			#say $ofas join "\n", ">".$id, $seq;
			$is_mariner = 1 if $seq =~ /ta/i;
		    }
		    unlink $tmp;
		    
		    #$is_mariner = 1 if $tsd =~ /ta/i;
		}
	    }
	    $mar_feats .= join "\t", @feats, "\n";
	}
	if ($is_mariner) {
	    chomp $mar_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out $mar_feats;
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdoms++ if $has_pdoms;
	}
	undef $mar_feats;
	$is_mariner = 0;
	$has_pdoms  = 0;
    }
    close $out;

    #if ($count > 0) {
	#say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = $stat->mean;
	my $count = $stat->count;
    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub find_hat {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $hat_feats;
    my $is_hat = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_hAT.gff3");
    open my $out, '>>', $outfile;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    $has_pdoms = 1 if $feats[2] =~ /protein_match/;
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		$is_hat = 1 if $tsd_len == 8;
	    }
	    $hat_feats .= join "\t", @feats, "\n";
	}
	if ($is_hat) {
	    chomp $hat_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out $hat_feats;
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdoms++ if $has_pdoms;
	}
	undef $hat_feats;
	$is_hat = 0;
	$has_pdoms = 0;
    }
    close $out;

    #if ($count > 0) {
    #say STDERR join "\t", "hat_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = $stat->mean;
	my $count = $stat->count;
    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub find_mutator {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $mut_feats;
    my $is_mutator = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_mutator.gff3");
    open my $out, '>>', $outfile;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    $has_pdoms = 1 if $feats[2] =~ /protein_match/;
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		if ($tsd_len >= 8 && $tsd_len <= 11) {
		    $is_mutator = 1;
		}
	    }
	    $mut_feats .= join "\t", @feats, "\n";
	}
	if ($is_mutator) {
	    chomp $mut_feats;
	    say $out $mut_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdoms++ if $has_pdoms;
	}
	undef $mut_feats;
	$is_mutator = 0;
	$has_pdoms = 0;
    }
    close $out;

    #if ($count > 0) {
	#say STDERR join "\t", "mutator_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = $stat->mean;
	my $count = $stat->count;
    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub find_cacta {
    my $self = shift;
    my ($feature, $header) = @_;
    my $fasta = $self->genome;
    my $gff   = $self->gff;
    
    my @lengths;
    my $cac_feats;
    my $is_cacta = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_cacta.gff3");
    open my $out, '>>', $outfile;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    $has_pdoms = 1 if $feats[2] =~ /protein_match/;
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		if ($tsd_len >= 2 && $tsd_len <= 3) {
		    #my $tir_elem = substr $hash->{$loc}, $s, $len;
		    my $tmp = $tir.".fasta";
		    my $cmd = "samtools faidx $fasta $region:$feats[3]-$feats[4] > $tmp";
		    $self->run_cmd($cmd);
		    my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
		    while (my $seqobj = $seqio->next_seq) {
			my $seq = $seqobj->seq;
			#$seq =~ s/.{60}\K/\n/g;
			#say $ofas join "\n", ">".$id, $seq;
			$is_cacta = 1 if $seq =~ /^cact[ag]/i;
		    }
		    unlink $tmp;
		    
		    #$is_cacta = 1 if $tir_elem =~ /^cact[ag]/i;
		}
	    }
	    $cac_feats .= join "\t", @feats, "\n";
	}
	if ($is_cacta) {
	    chomp $cac_feats;
	    say $out $cac_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdoms++ if $has_pdoms;
	}
	undef $cac_feats;
	$is_cacta = 0;
	$has_pdoms = 0;
    }
    close $out;

    #if ($count > 0) {
	#say STDERR join "\t", "cacta_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = $stat->mean;
	my $count = $stat->count;
    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub write_unclassified_tirs {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $unc_feats;
    my $is_unclass = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_unclassified.gff3");
    open my $out, '>>', $outfile;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+/=/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/=$//;
	    $feats[8] =~ s/=\;/;/g;
	    $feats[8] =~ s/\"//g;
	    $has_pdoms = 1 if $feats[2] =~ /protein_match/;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out join "\t", @feats;
	}
	delete $feature->{$tir};
	push @lengths, $len;
	$pdoms++ if $has_pdoms;
	$has_pdoms = 0;
    }
    close $out;

    #if ($count > 0) {
	#say STDERR join "\t", "unclassified_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	my $stat = Statistics::Descriptive::Full->new;
	$stat->add_data(@lengths);
	my $min   = $stat->min;
	my $max   = $stat->max;
	my $mean  = $stat->mean;
	my $count = $stat->count;
    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

__PACKAGE__->meta->make_immutable;

1;
