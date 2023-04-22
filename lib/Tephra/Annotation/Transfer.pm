package Tephra::Annotation::Transfer;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use File::Find;
use File::Basename;
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use List::Util          qw(sum);
use Carp                qw(croak);
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Annotation::Transfer - Transfer annotations from a reference set to Tephra annotations

=head1 VERSION

Version 0.14.0

=cut

our $VERSION = '0.14.0';
$VERSION = eval $VERSION;

#
# methods
#
sub transfer_annotations {
    my $self = shift;

    my $blast_report = $self->process_blast_args;
    my ($combined, $id_map) = $self->parse_blast($blast_report);
    $self->write_annotations($combined, $id_map);
    unlink $blast_report;
}

sub parse_blast {
    my $self = shift;
    my $blast_hpid = $self->blast_hit_pid;
    my $blast_hcov = $self->blast_hit_cov;
    my $blast_hlen = $self->blast_hit_len;
    my ($blast_report) = @_;

    my $perc_cov = sprintf("%.2f",$blast_hcov/100);

    my %matches;
    open my $in, '<', $blast_report or die "\n[ERROR]: Could not open file: $blast_report\n";
    while (my $line = <$in>) {
	chomp $line;
	my ($queryid, $hitid, $pid, $hitlen, $mmatchct, $gapct, 
	    $qstart, $qend, $hstart, $hend, $evalue, $score) = split /\t/, $line;
	my ($start, $end) = ($queryid =~ /(\d+)-?_?(\d+)$/);
	my $qlen = $end - $start + 1;
	if ($hitlen >= $blast_hlen && $hitlen >= ($qlen * $perc_cov) && $pid >= $blast_hpid) {
	    push @{$matches{$queryid}{$score}}, 
	        join "||", $queryid, $hitid, $pid, $hitlen, $mmatchct, $gapct, $qstart, $qend, $hstart, $hend, $evalue, $score;
	}
    }
    close $in;

    my (%hits, %id_map);
    for my $query (keys %matches) {
	my $family = ($query =~ /(^\w+_family\d+)_/) ? $1 : $query;
	push @{$id_map{ $family }}, $query;
	for my $score (reverse sort { $a <=> $b } keys %{$matches{$query}}) {
	    my $hitct = @{$matches{$query}{$score}};
	    for my $match (@{$matches{$query}{$score}}) {
		my @f = split /\|\|/, $match;
		my $hit_fam = $f[1];
		$hit_fam =~ s/I$|_I$|LTR$|_LTR$//; # RepBase
		$hit_fam =~ s/_AC.*//;             # maizeTEdb
		$hit_fam =~ s/-/_/g;               # sunflowerTEdb
		$hit_fam =~ s/\d+$//;              # sunflowerTEdb
		$hit_fam =~ s/(\w{3}_\w+?)_.*/$1/; # sunflowerTEdb
		push @{$hits{$family}{$hit_fam}}, join "||", @f[0, 1, 2, 3, 11];
	    }
	}
    }

    #dd \%hits and exit;
    my (%fam_map, %elements, @hits, @matches);
    for my $family (sort keys %hits) {
	for my $hit_type (keys %{$hits{$family}}) {
	    my $hit_ct = @{$hits{$family}{$hit_type}};
	    for my $match (@{$hits{$family}{$hit_type}}) {
		my @matches = split /\|\|/, $match;
		#RLC_family0_LTR_retrotransposon63617_9_16309492_16318958||RLC_ji_AC210731-10832||96.04||7645||12414
		push @hits, $matches[4];
	    }
	    my $score_sum = sum(@hits);
	    $fam_map{$family}{$hit_type} = { hit_ct => $hit_ct, score_sum => $score_sum };
	    @hits = ();
	    %elements = ();
	}
    }

    my $ct = 0;
    my ($max_ct, $max_score, $best_hit);
    my (%combined, %rev_fam_map, @all_scores, @all_hitcts);
    for my $family (sort keys %fam_map) { 
	my $fam_hits = keys %{$fam_map{$family}};
	for my $hit (keys %{$fam_map{$family}}) {
	    my ($hit_sum, $scores) = @{$fam_map{$family}{$hit}}{qw(hit_ct score_sum)};
	    if ($ct > 0) {
		if ($hit_sum > $max_ct && $scores > $max_score) {
		    $max_ct    = $hit_sum;
		    $max_score = $scores;
		    $best_hit  = $hit;
		}
	    }
	    else {
		$ct++;
		$max_ct    = $hit_sum;
		$max_score = $scores;
		$best_hit  = $hit;
	    }
	}
	push @{$combined{$best_hit}}, $family;
	$best_hit  = '';
	$max_score = 0;
	$max_ct    = 0;
	
    }

    return (\%combined, \%id_map);
}

sub write_annotations {
    my $self = shift;
    my $fasta   = $self->infile;
    my $outfile = $self->outfile;
    my ($combined, $id_map) = @_;

    my $index = $self->index_ref($fasta);
    my $kseq = Bio::DB::HTS::Kseq->new($fasta);
    my $iter = $kseq->iterator;

    open my $outfh, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";

    #{
	#"DTC_ZM00004_consensus" => ["RLX_singleton_family2330"],
  
    my ($total, $transferred, $not_transferred) = (0, 0, 0);
    my %seen;
    for my $family_match (keys %$combined) {
	for my $element (@{$combined->{$family_match}}) {
	    my $family = ($element =~ /(^\w+_family\d+)_/) ? $1 : $element;
	    my $id_arr = $id_map->{$family};
	    for my $elemid (@$id_arr) {
		$transferred++;
		$seen{ $elemid } = 1;
		my ($seq, $length) = $index->get_sequence($elemid);
		$seq =~ s/.{60}\K/\n/g;
		$elemid =~ s/${family}_//;
		say $outfh join "\n", ">".$family_match."_".$elemid, $seq;
	    }
	}
    }

    while (my $seqio = $iter->next_seq) {
	$total++;
	my $id = $seqio->name;
	my $seq = $seqio->seq;
	$seq =~ s/.{60}\K/\n/g;
	unless (exists $seen{ $id }) {
	    $not_transferred++;
	    say $outfh join "\n", ">".$id, $seq;
	}
    }
    close $outfh;

    
    my $perc_transferred = sprintf("%.2f", $transferred/$total * 100);
    
    say STDERR "=====> $perc_transferred% ($transferred/$total) annotations transferred from reference.";
    #say STDERR "=====> Total: $total";
    #say STDERR "=====> Transferred: $transferred";
    #say STDERR "=====> Not transferred: $not_transferred";
    return;
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

    perldoc Tephra::Annotation::Transfer


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
