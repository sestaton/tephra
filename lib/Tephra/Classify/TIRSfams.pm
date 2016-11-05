package Tephra::Classify::TIRSfams;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::GFF3::LowLevel qw(gff3_format_feature);
use IPC::System::Simple qw(capture);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use Path::Class::File;
use Try::Tiny;
use Cwd;
use Carp 'croak';
use Tephra::Config::Exe;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

=head1 NAME

Tephra::Classify::TIRSams - Classify TIR transposons into superfamilies

=head1 VERSION

Version 0.04.2

=cut

our $VERSION = '0.04.2';
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
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_tc1-mariner.gff3');
    my $fas     = File::Spec->catfile($path, $name.'_tc1-mariner.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_tc1-mariner_domain_org.tsv');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\nERROR: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";;
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
	my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
	for my $tir_feature (@{$feature->{$rep_region}}) {
	    if ($tir_feature->{type} eq 'protein_match') {
		$has_pdoms = 1;
		my $pdom_name = $tir_feature->{attributes}{name}[0];
		push @all_pdoms, $pdom_name;
	    }

	    if ($tir_feature->{type} eq 'target_site_duplication') {
		my ($seq_id, $tsd_start, $tsd_end) = @{$tir_feature}{qw(seq_id start end)};
		my $tsd_len = $tsd_end - $tsd_start + 1;
		if ($tsd_len == 2) {
		    my $seq = $self->subseq($index, $seq_id, undef, $tsd_start, $tsd_end, undef);
		    $is_mariner = 1 if $seq =~ /ta/i;
		}
	    }

	    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		my $elem_id = $tir_feature->{attributes}{ID}[0];
		($seq_id, $source, $start, $end, $strand) 
		    = @{$tir_feature}{qw(seq_id source start end strand)};
		my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		my $id = join "_", 'DTT', $elem_id, $seq_id, $start, $end;

		$lines .= join "\n", ">".$id, $seq;
		$len = $end - $start + 1;
            }
	    my $gff3_str = gff3_format_feature($tir_feature);
	    $mar_feats .= $gff3_str;
	}

	if ($is_mariner) {
	    chomp $mar_feats;
	    say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
	    say $out $mar_feats;
	    delete $feature->{$rep_region};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;

	    say $faout $lines;
	}	    
	undef $mar_feats;
	undef $pdom_org;
	undef $lines;

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
	unlink $outfile, $fas;	
    }
}

sub find_hat {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff   = $self->gff;
    my $fasta = $self->genome;
    
    my @lengths;
    my $hat_feats;
    my $is_hat = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_hAT.gff3');
    my $fas     = File::Spec->catfile($path, $name.'_hAT.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_hAT_domain_org.tsv');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\nERROR: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $id, $seq_id, $source, $start, $end, $strand);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
            if ($tir_feature->{type} eq 'protein_match') {
                $has_pdoms = 1;
                my $pdom_name = $tir_feature->{attributes}{name}[0];
                push @all_pdoms, $pdom_name;
            }

            if ($tir_feature->{type} eq 'target_site_duplication') {
                my ($seq_id, $tsd_start, $tsd_end) = @{$tir_feature}{qw(seq_id start end)};
                my $tsd_len = $tsd_end - $tsd_start + 1;
		$is_hat = 1 if $tsd_len == 8;
            }

	    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
		my $elem_id = $tir_feature->{attributes}{ID}[0];
                ($seq_id, $source, $start, $end, $strand) 
                    = @{$tir_feature}{qw(seq_id source start end strand)};
                my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
                $id = join "_", 'DTA', $elem_id, $seq_id, $start, $end;

                $lines .= join "\n", ">".$id, $seq;
                $len = $end - $start + 1;
            }
            my $gff3_str = gff3_format_feature($tir_feature);
            $hat_feats .= $gff3_str;
	}

	if ($is_hat) {
	    chomp $hat_feats;
	    say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
            say $out $hat_feats;
            delete $feature->{$rep_region};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	    
	    say $faout $lines;
	}
	undef $hat_feats;
	undef $pdom_org;
	undef $lines;

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

    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "hat_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile, $fas;
    }
}

sub find_mutator {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff   = $self->gff;
    my $fasta = $self->genome;

    my @lengths;
    my $mut_feats;
    my $is_mutator = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_mutator.gff3');
    my $fas     = File::Spec->catfile($path, $name.'_mutator.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_mutator_domain_org.tsv');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\nERROR: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $id, $seq_id, $source, $start, $end, $strand);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
            if ($tir_feature->{type} eq 'protein_match') {
                $has_pdoms = 1;
                my $pdom_name = $tir_feature->{attributes}{name}[0];
                push @all_pdoms, $pdom_name;
            }

            if ($tir_feature->{type} eq 'target_site_duplication') {
                my ($seq_id, $tsd_start, $tsd_end) = @{$tir_feature}{qw(seq_id start end)};
                my $tsd_len = $tsd_end - $tsd_start + 1;
		if ($tsd_len >= 8 && $tsd_len <= 11) {
		    $is_mutator = 1;
		}
	    }

	    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
                my $elem_id = $tir_feature->{attributes}{ID}[0];
                ($seq_id, $source, $start, $end, $strand) 
                    = @{$tir_feature}{qw(seq_id source start end strand)};
                my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
                my $id = join "_", 'DTM', $elem_id, $seq_id, $start, $end;

                $lines .= join "\n", ">".$id, $seq;
                $len = $end - $start + 1;
            }
            my $gff3_str = gff3_format_feature($tir_feature);
            $mut_feats .= $gff3_str;
        }

        if ($is_mutator) {
            chomp $mut_feats;
            say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
            say $out $mut_feats;
            delete $feature->{$rep_region};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;

	    say $faout $lines;
	}
	undef $mut_feats;
	undef $pdom_org;
	undef $lines;

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
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "mutator_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";	
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile, $fas;
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
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_cacta.gff3');
    my $fas     = File::Spec->catfile($path, $name.'_cacta.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_cacta_domain_org.tsv');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\nERROR: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand, $tir_len);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
            if ($tir_feature->{type} eq 'protein_match') {
                $has_pdoms = 1;
                my $pdom_name = $tir_feature->{attributes}{name}[0];
                push @all_pdoms, $pdom_name;
            }

	    if ($tir_feature->{type} eq 'terminal_inverted_repeat') {
		my ($tir_start, $tir_end) = @{$tir_feature}{qw(start end)};                            
		$tir_len = $tir_end - $tir_start + 1; 
	    }

	    if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
                my $elem_id = $tir_feature->{attributes}{ID}[0];
                ($seq_id, $source, $start, $end, $strand) 
                    = @{$tir_feature}{qw(seq_id source start end strand)};
                my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
		if ($seq =~ /^cact(?:a|g)?|cact(?:a|g)?$/i) {
		    # Lewin, 1997 http://www.plantphysiol.org/content/132/1/52.full
		    # provides this TIR length definition, but it seems to remove all predictions
		    $is_cacta = 1; # if $tir_len >= 10 && $tir_len <= 28;
		}
		
                my $id = join "_", 'DTC', $elem_id, $seq_id, $start, $end;

                $lines .= join "\n", ">".$id, $seq;
                $len = $end - $start + 1;
            }
            my $gff3_str = gff3_format_feature($tir_feature);
            $cac_feats .= $gff3_str;
        }

        if ($is_cacta) {
            chomp $cac_feats;
            say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
            say $out $cac_feats;
            delete $feature->{$rep_region};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;

	    say $faout $lines;
	}
	undef $cac_feats;
	undef $pdom_org;
	undef $lines;

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
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "cacta_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile, $fas;
    }
}

sub write_unclassified_tirs {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff   = $self->gff;
    my $fasta = $self->genome;

    my @lengths;
    my $unc_feats;
    my $is_unclass = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_unclassified.gff3');
    my $fas     = File::Spec->catfile($path, $name.'_unclassified.fasta');
    my $domoutfile = File::Spec->catfile($path, $name.'_unclassified_domain_org.tsv');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\nERROR: Could not open file: $fas\n";
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";
    say $out $header;

    my ($len, $lines, $seq_id, $source, $start, $end, $strand);
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$feature) {
        my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $tir_feature (@{$feature->{$rep_region}}) {
            if ($tir_feature->{type} eq 'protein_match') {
                $has_pdoms = 1;
                my $pdom_name = $tir_feature->{attributes}{name}[0];
                push @all_pdoms, $pdom_name;
            }

            if ($tir_feature->{type} eq 'terminal_inverted_repeat_element') {
                my $elem_id = $tir_feature->{attributes}{ID}[0];
                ($seq_id, $source, $start, $end, $strand) 
                    = @{$tir_feature}{qw(seq_id source start end strand)};
                my $seq = $self->subseq($index, $seq_id, $elem_id, $start, $end, undef);
                my $id = join "_", 'DTX', $elem_id, $seq_id, $start, $end;

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
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "unclassified_tir_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile, $fas;
    }
}

sub subseq {
    my $self = shift;
    my ($index, $loc, $elem, $start, $end, $out) = @_;

    my $location = "$loc:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);
    croak "\nERROR: Something went wrong. This is a bug, please report it.\n"
        unless $length;

    $seq =~ s/.{60}\K/\n/g;

    return $seq;
}

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

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
