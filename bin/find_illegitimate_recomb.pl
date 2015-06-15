#!/usr/bin/env perl

# TODO: - Write flanking gap seqs to separate file for each input seq in alignment
#       - Take search space as option as well as PID for flanking gap  sequence matches
#       - Printing the pairwise alignments is still messed up so we are just going to bypass that routine for now
#       - Remove use of BioPerl; add global alignment methods

use 5.010;
use strict;
use warnings;
use Cwd;
use Bio::Seq;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::Tools::Run::StandAloneBlast;
use Statistics::Descriptive;
use Getopt::Long;
use File::Basename;
use File::Temp;
#use Text::LevenshteinXS;       # This is exactly what I want: For string comparison to get the edit distance (num deletions, insertions)
#use String::Approx;            # This may also be useful: For string matching
#use Data::Dumper

my $aln_file;
my $outfile;
my $statsfile;
my $dr_pid;
my $gap_stats;
my $help;

my $pos = 0;
my $gap = 0;
my $del = 0;
my @indels;
my @flanking_seqs;

GetOptions(
	   'i|infile=s'     => \$aln_file,
	   'o|outfile=s'    => \$outfile,
	   's|statsfile=s'  => \$statsfile,
	   'p|repeat_pid=i' => \$dr_pid,
           'h|help'         => \$help,
	   );

usage() and exit(0) if $help;

usage() and exit(1) if !$aln_file;

$dr_pid = defined($dr_pid) ? $dr_pid : '10';
my $cwd = getcwd();
my $statstmp = $statsfile.".tmp";
open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
open my $stats_out_tmp, '>>', $statstmp or die "\nERROR: Could not open file: $statstmp\n";
say $stats_out_tmp "#LTR_retro_name\tTotal_num_gap_sites\tNum_diff_gap_sizes\tPercent_gap_sites\tMean_gap_size\tMin_gap_size\tMax_gap_size";

my ($seqs_in_aln,$count) = split_aln($aln_file);

for my $fas (@$seqs_in_aln) {
    my $aln_in = Bio::AlignIO->new(-fh => \*$fas, -format => 'fasta');

    my ($fname, $fpath, $fsuffix) = fileparse($fas, qr/\.[^.]*/);
    my $seq_out = $fname;
    $seq_out .= "_gap_flanking_sequences.fasta";
    open my $each_out, '>>', $seq_out or die "\nERROR: Could not open file: $seq_out\n";

    while ( my $aln = $aln_in->next_aln() ) {
	my $aln_len = $aln->length;
	for my $hash (@{$aln->gap_col_matrix}) {
	    for my $key (keys %$hash) {
		$pos++; 
		if ($hash->{$key} eq '1') { # gap columns are coded as '1'
		    $gap++;
		    push @indels, $pos; 
		}   
	    }
	}
	
	my ($indel_lengths, $indel_ranges) = get_indel_range(@indels);   
	@indels = ();

	$gap_stats = get_stats($fas,$pos,$gap,$indel_lengths,$fname) if $statsfile;
	
    	for my $line (@$indel_ranges) {
	    $del++;
	    my ($indel_len, $indel_spos, $indel_epos) = split(/\t/,$line); 

	    next unless $indel_len >= 10; # Only analyze repeats flanking deletions > 10 bp; Ma et al. 2004, Genome Research
	    
	    # What we want to do is search 20 bp on either side of the gap. 
	    # This is based on the 1-15 direct repeats found in Arabidopsis... (Devos et al. 2002)
	    # Sunflower may be different, so we search 20 bp.
	    my $upstream_spos   = $indel_spos - 20;
	    my $upstream_epos   = $indel_spos - 1;
	    my $downstream_spos = $indel_epos + 1;
	    my $downstream_epos = $indel_epos + 20;
        
	    if ($upstream_spos < 1 || $upstream_epos > $aln_len) {                                                    
		say "Deletion $del has a flanking repeat out of bounds.";
		next;
	    }
	    if ($downstream_epos > $aln_len) {
		say "Deletion $del has the downstream flanking repeat out of bounds.";
		last;
	    }
	    
	    for my $seq ($aln->each_seq) { # $seq "isa" a Bio::LocatableSeq object (which is a part of Bio::PrimarySeq)
		my $upstr_seq   = $seq->subseq($upstream_spos,   $upstream_epos);
		my $downstr_seq = $seq->subseq($downstream_spos, $downstream_epos);
		
		my $upstr_id   = $seq->id."_upstr-del-".$del."_".$upstream_spos."-".$upstream_epos;
		my $downstr_id = $seq->id."_downstr-del-".$del."_".$downstream_spos."-".$downstream_epos;
		
		my $upstream_seqobj = Bio::Seq->new(-seq      => $upstr_seq,
						    -id       => $upstr_id,
						    -alphabet => 'dna');
		
		my $downstream_seqobj = Bio::Seq->new(-seq      => $downstr_seq,
						      -id       => $downstr_id,
						      -alphabet =>'dna');
		
		# Here is where we do the comparison
		blast_compare($indel_len, $upstream_seqobj, $downstream_seqobj, $cwd, $each_out);
		
		#//////// could call a sub or just do the Levenshtein distance($upstr_seq, $downstr_seq) here \\\\\\\\\
	    }       	    
	}
    }
    $pos = 0;
    $del = 0;
    close $each_out;
    collate $seq_out,$out;
    collate $gap_stats,$stats_out_tmp;
    unlink $fas;
    unlink $seq_out;
    unlink $gap_stats;
}

close $out;
close $stats_out_tmp;

collate_gap_stats($statstmp,$statsfile);     
unlink $statstmp;

exit;
#
# methods
#
sub split_aln {
    my ($input) = @_;

    my ($iname, $ipath, $isuffix) = fileparse($input, qr/\.[^.]*/);

    my $seq_in  = Bio::SeqIO->new(-file  => $input,
                                  -format => 'fasta');
    my $count = 0;
    my $fcount = 1;
    my @split_files;
    $iname =~ s/\.fa.*//;     # clean up file name like seqs.fasta.1
    
    my $tmpiname = $iname."_".$fcount."_XXXX";
    my $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                 DIR      => $cwd,
                                 UNLINK   => 0,
				 SUFFIX   => ".fasta");
    
    my $seq_out = Bio::SeqIO->new(-file   => ">$fname", 
                                  -format => 'fasta');

    push @split_files, $fname;
    while (my $seq = $seq_in->next_seq) {
        if ($count % 1 == 0 && $count > 0) {
            $fcount++;
            $tmpiname = $iname."_".$fcount."_XXXX";
            $fname = File::Temp->new( TEMPLATE => $tmpiname,
                                      DIR      => $cwd,
                                      UNLINK   => 0,
				      SUFFIX   => ".fasta");

            $seq_out = Bio::SeqIO->new(-file   => ">$fname", 
                                       -format => 'fasta');

            push(@split_files,$fname);
        }
        $seq_out->write_seq($seq);
        $count++;
    }

    return (\@split_files, $count);
}

sub get_indel_range {
    my @indels = @_;
    my @indel_ranges;
    my @indel_lengths;

    # The algorithm below is based on 
    # something I found on stackoverflow
    # for converting an array with sequential 
    # numbers into an array with ranges
    # http://stackoverflow.com/questions/117691
    my $gap_range = [ $indels[0] ];

    for my $indel (@indels[1..$#indels]) {
        if ($indel == $gap_range->[-1] + 1) {
            push @$gap_range, $indel;
        }
        else {
            my $gap_length = ($gap_range->[-1] - $gap_range->[0]) + 1;
            push @indel_lengths, $gap_length;
            push @indel_ranges, $gap_length."\t".$gap_range->[-1]."\t".$gap_range->[0]."\n";
            $gap_range = [ $indel ];
        }
    }
    return (\@indel_lengths, \@indel_ranges);
}

sub blast_compare {
    my ($indel_len, $upstream_seqobj, $downstream_seqobj, $cwd, $each_out) = @_;
	
    my $bl2_out = File::Temp->new(TEMPLATE => "bl2seq_out_XXXX", DIR => $cwd);
    my $bl2_out_fname = $bl2_out->filename;

    my @params = (-F => 'F', -p => 'blastn', -W => 4, -o => $bl2_out_fname);
    my $factory = Bio::Tools::Run::StandAloneBlast->new(@params);
								    		
    $factory->bl2seq($upstream_seqobj, $downstream_seqobj);
    
    my $bl2seq_report = Bio::SearchIO->new(-file => $bl2_out_fname, -format => 'blast');
    
    while ( my $result = $bl2seq_report->next_result ) {
	my $query      = $result->query_name();
	my $qlen       = $result->query_length();
	
	while ( my $hit = $result->next_hit ) {
	    my $hitid    = $hit->name();
	    
	    while ( my $hsp = $hit->next_hsp ) {		
		my $hsplen    = $hsp->length('total');
		my $hstart    = $hsp->start('hit');
		my $hstop     = $hsp->end('hit'); 
		my $qstart    = $hsp->start('query');
		my $qstop     = $hsp->end('query');
		my $qstring   = $hsp->query_string;
		my $hstring   = $hsp->hit_string;
		my $homstring = $hsp->homology_string; 
		my $hpid      = $hsp->percent_identity;
		
		if ( $hsplen > 2 && $hpid >= $dr_pid ) {    # $dr_pid ("direct repeat % identity") is 10 or 2/20
		    my $qustr = $upstream_seqobj->seq;
                    my $hitstr = $downstream_seqobj->seq;
		
		    #if ($hstring !~ /$downstream_seqobj->seq/i) {
			#print "Why the lack of match?\n";
			#next;
		    #}
		    my $match_hit_string = uc($hstring);
		
		    if ($hitstr !~ m/$match_hit_string/i) {
			#print "$match_hit_string\t$hitstr\n"; #BioPerl bug? Why is $qstring being reported for $hstring sometimes?
			next;                                  #This just ensures the $hstring came from the actual Hit sequence
		    }
		    say $each_out join "\t", "Query_ID","Hit_ID","HSP_len","Hit_start","Hit_stop","Query_start","Query_stop","HSP_PID";
		    say $each_out join "\t", $query,$hitid,$hsplen,$hstart,$hstop,$qstart,$qstop,$hpid,"\n";     
	

		    say $each_out "Gap length        : $indel_len";
		    say $each_out "Query len         : $qlen";
		    say $each_out "Hit len           : ",$hit->length;
		    say $each_out "Query string      : $qustr\n";
		    say $each_out "Query match string: ",uc($qstring);
		    say $each_out "Homology string   : $homstring";
		    say $each_out "Hit match string  : ",uc($hstring),"\n";
		    say $each_out "Hit string        : $hitstr\n";
		    
		}
	    }
	}
    }
}

sub get_stats {
    my ($fas, $pos, $gap, $indel_lengths, $fname) = @_;

    my $gap_stats = $fname."_gap_stats.txt";
    open my $statsout, '>', $gap_stats or die "\nERROR: Could not open file: $gap_stats\n";
    my $stat = Statistics::Descriptive::Full->new;

    $stat->add_data(@$indel_lengths);

    my ($fasname, $faspath, $fassuffix) = fileparse($fas, qr/\.[^.]*/);
    my $count = $stat->count;
    my $gap_percent = sprintf("%.2f",$gap/$pos);
    my $mean = sprintf("%.2f",$stat->mean);
    my $min = $stat->min;
    my $max = $stat->max;

    say $statsout "$fasname\t$gap\t$count\t$gap_percent\t$mean\t$min\t$max";
    
    close $statsout;
    return $gap_stats;
}

sub collate {
    my ($file_in, $fh_out) = @_;
    open my $fh_in, '<', $file_in or die "\nERROR: Could not open file: $file_in\n";
    while (<$fh_in>) {
	print $fh_out $_;
    }
    close $fh_in;
}

sub collate_gap_stats {
    my ($statstmp, $statsfile) = @_;

    open my $gap_stats_fh_in, '<', $statstmp or die "\nERROR: Could not open file: $statstmp\n";
    open my $gap_stats_fh_out, '>', $statsfile or die "\nERROR: Could not open file: $statsfile\n";

    my (@repeat_names, @total_gap_char, @diff_gap_sizes, @gap_char_perc, @mean_gap_size, @min_gap_size, @max_gap_size);

    while (<$gap_stats_fh_in>) {
	if (/^\#/) {
	    print $gap_stats_fh_out $_;
	} 
	else {
	    my @all_gap_stats = split /\t/, $_;
	    push @repeat_names,   $all_gap_stats[0];
	    push @total_gap_char, $all_gap_stats[1];
	    push @diff_gap_sizes, $all_gap_stats[2];
	    push @gap_char_perc,  $all_gap_stats[3];
	    push @mean_gap_size,  $all_gap_stats[4];
	    push @min_gap_size,   $all_gap_stats[5];
	    push @max_gap_size,   $all_gap_stats[6];
 
	    print $gap_stats_fh_out $_;
	}
    }

    my $fam_name = pop @repeat_names;
    $fam_name =~ s/\_\d+\_.*//;
    $fam_name =~ s/^all\_//;

    my $total_gap_char_stats = Statistics::Descriptive::Full->new;
    my $gap_size_stats       = Statistics::Descriptive::Full->new;
    my $gap_char_perc_stats  = Statistics::Descriptive::Full->new;
    my $mean_gap_size_stats  = Statistics::Descriptive::Full->new;
    my $min_gap_size_stats   = Statistics::Descriptive::Full->new;
    my $max_gap_size_stats   = Statistics::Descriptive::Full->new;

    $total_gap_char_stats->add_data(@total_gap_char);
    $gap_size_stats->add_data(@diff_gap_sizes);
    $gap_char_perc_stats->add_data(@gap_char_perc);
    $mean_gap_size_stats->add_data(@mean_gap_size);
    $min_gap_size_stats->add_data(@min_gap_size);
    $max_gap_size_stats->add_data(@max_gap_size);

    my $grand_mean_fam_count = $total_gap_char_stats->count;
    my $grand_gap_char_mean  = sprintf("%.2f",$total_gap_char_stats->mean);
    my $grand_gap_size_mean  = sprintf("%.2f",$gap_size_stats->mean);
    my $grand_gap_char_perc  = sprintf("%.2f",$gap_char_perc_stats->mean);
    my $grand_mean_gap_size  = sprintf("%.2f",$mean_gap_size_stats->mean);
    my $grand_gap_size_min   = sprintf("%.2f",$min_gap_size_stats->mean);
    my $grand_gap_size_max   = sprintf("%.2f",$max_gap_size_stats->mean);

    my $grand_gap_char_sd      = sprintf("%.2f",$total_gap_char_stats->standard_deviation);
    my $grand_gap_size_sd      = sprintf("%.2f",$gap_size_stats->standard_deviation);
    my $grand_gap_char_perc_sd = sprintf("%.2f",$gap_char_perc_stats->standard_deviation);
    my $grand_gap_size_mean_sd = sprintf("%.2f",$mean_gap_size_stats->standard_deviation);
    my $grand_gap_size_min_sd  = sprintf("%.2f",$min_gap_size_stats->standard_deviation);
    my $grand_gap_size_max_sd  = sprintf("%.2f",$max_gap_size_stats->standard_deviation);

    say $gap_stats_fh_out "\nFamily_name\tTotal_fam_gap_count\tMean_fam_gap_count ".
                          "(stddev)\tMean_fam_gap_size (stddev)\tMean_fam_gap_perc ".
	                  "(stdev)\tMean_fam_gap_size (stddev)\tMean_gap_min_size ".
                          "(stddev)\tMean_gap_max_size (stddev)";

    say $gap_stats_fh_out $fam_name."\t".$grand_mean_fam_count."\t".
	                  $grand_gap_char_mean." (".$grand_gap_char_sd.")"."\t".
			  $grand_gap_size_mean." (".$grand_gap_size_sd.")"."\t".
			  $grand_gap_char_perc." (".$grand_gap_char_perc_sd.")"."\t".
			  $grand_mean_gap_size." (".$grand_gap_size_mean_sd.")"."\t".
			  $grand_gap_size_min." (".$grand_gap_size_min_sd.")"."\t".
			  $grand_gap_size_max." (".$grand_gap_size_max_sd.")"."\n";


    close $gap_stats_fh_in;
    close $gap_stats_fh_out;
}

sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i seqs.fas -o blast_result 

Required:
    -i|infile        :    An alignment file in fasta format.
    -o|outfile       :    File name to write the extracted sequences to.

Options:
    -s|statsfile     :    A file to write alignment stats to.
    -p|repeat_pid    :    The percent identity threshold for retaining repeats that flank gaps. (Default: 10, which is a 2 bp match).

END
}
