#!/usr/bin/env perl

## Take the gff output file from gt tirvish and classify TEs into superfamilies.
##
## NB: It is important to first filter out retroelements.

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Statistics::Descriptive;
use Bio::Tools::GFF;
use List::UtilsBy qw(nsort_by);
use Data::Dump;
use experimental 'signatures';

my $usage = "$0 gff fasta";
my $gff   = shift or die $usage;
my $fasta = shift or die $usage;
my $header;

my $hash = seq_to_hash($fasta);

open my $in, '<', $gff;
while (<$in>) {
    chomp;
    if (/^#/) {
	$header .= $_."\n";
    }
    else {
	last;
    }
}
close $in;
chomp $header;

my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

my ($start, $end, $region, %feature);
while (my $feature = $gffio->next_feature()) {
    if ($feature->primary_tag eq 'repeat_region') {
	my @string = split /\t/, $feature->gff_string;
	($region) = ($string[8] =~ /ID=?\s+?(repeat_region\d+)/);
	($start, $end) = ($feature->start, $feature->end);
    }
    next $feature unless defined $start && defined $end;
    if ($feature->primary_tag ne 'repeat_region') {
	if ($feature->start >= $start && $feature->end <= $end) {
	    push @{$feature{$region.".".$start.".".$end}}, join "||", split /\t/, $feature->gff_string;
	}
    }
}

my $all_ct = (keys %feature);
find_tc1_mariner(\%feature, $header, $hash, $gff);
my $tc1_ct = (keys %feature);
find_hat(\%feature, $header, $hash, $gff);
my $hat_ct = (keys %feature);
find_mutator(\%feature, $header, $hash, $gff);
my $mut_ct = (keys %feature);
find_cacta(\%feature, $header, $hash, $gff);
my $cacta_ct = (keys %feature);
write_unclassified_tirs(\%feature, $header, $hash, $gff);
my $rem_ct = (keys %feature);

say STDERR join "\t", "all", "after_tc1", "after_hat", "after_mut", "after_cacta", "after_rem";
say STDERR join "\t", $all_ct, $tc1_ct, $hat_ct, $mut_ct, $cacta_ct, $rem_ct;
#
# methods
#
sub find_tc1_mariner ($feature, $header, $hash, $gff) {
    my @lengths;
    my $mar_feats;
    my $is_mariner = 0;
    my $has_pdoms  = 0;
    my $pdoms      = 0;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_tc1-mariner.gff3";
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
		    my $tsd = substr $hash->{$loc}, $feats[3], $tsd_len;
		    $is_mariner = 1 if $tsd =~ /ta/i;
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

    say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub find_hat ($feature, $header, $hash, $gff) {
    my @lengths;
    my $hat_feats;
    my $is_hat = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_hAT.gff3";
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
    
    say STDERR join "\t", "hat_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub find_mutator ($feature, $header, $hash, $gff) {
    my @lengths;
    my $mut_feats;
    my $is_mutator = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_mutator.gff3";
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

    say STDERR join "\t", "mutator_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub find_cacta ($feature, $header, $hash, $gff) {
    my @lengths;
    my $cac_feats;
    my $is_cacta = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_cacta.gff3";
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
                    my $tir_elem = substr $hash->{$loc}, $s, $len;
		    $is_cacta = 1 if $tir_elem =~ /^cact[ag]/i;
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

    say STDERR join "\t", "cacta_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub write_unclassified_tirs ($feature, $header, $hash, $gff) {
    my @lengths;
    my $unc_feats;
    my $is_unclass = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = $name."_unclassified.gff3";
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

    say STDERR join "\t", "unclassified_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
}

sub get_source {
    my ($ref) = @_;

    for my $feat (@$ref) {
	for my $rfeat (@$feat) {
	    my @feats = split /\|\|/, $rfeat;
	    return ($feats[0], $feats[1]);
	}
    }
}

sub seq_to_hash {
    my ($file) = @_;

    open my $in, '<', $file;
    my %hash;
    my @aux = undef;
    my ($name, $comm, $seq, $qual);

    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
        $hash{$name} = $seq;
    }
    close $in;

    return \%hash;
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}
