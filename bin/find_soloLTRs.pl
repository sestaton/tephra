#!/usr/bin/env perl

## TODO: - Create GFF of solo LTRs
##       - If a directory of models is given, just search with the models
##       - Improve file and path handling
##       - Add better command execution
##       - Consider removing shell commands

## USAGE: Just create exemplar for each fam and use that as input,
##        otherwise the number of files will explode with this script

## NB: This uses HMMER 2.3.2 because  HMMER3 does not currently support globabl alignments.
##
use 5.010;
use strict;
use warnings;
use Cwd;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename;
use File::Copy;
use File::Spec;
use Bio::AlignIO;
use Bio::SearchIO;

#
# define vars
#
my $in_dir; 
my $fas_dir;
my $hmm_path;
my $percentID_threshold;
my $quiet;
my $report;
my $global;
my $multi;
my $seq;
my $clean;
my $header;
my $calibrate;

GetOptions(# Required arguments
	   'i|indir=s'         => \$in_dir, 
           'f|fasdir=s'        => \$fas_dir,
	   'hmm_path=s'        => \$hmm_path,
	   'p|percentident=f'  => \$percentID_threshold,
 	   'quiet'             => \$quiet,
           'report'            => \$report,
	   'global'            => \$global,
	   'multi'             => \$multi,
	   'seq'               => \$seq,
	   'clean'             => \$clean,
	   'header'            => \$header,
	   'calibrate'         => \$calibrate,
	  );

#
# Check input
# 
if (!$in_dir || !$fas_dir) {
    print "\nERROR: Command line not parsed correctly. Exiting.\n\n";
    usage();
    exit(0);
} 
if (!$percentID_threshold) {
    warn "\nUsing 0.39 percent identity threshold between sequence and model matches because a value was not specified.\n" unless $quiet;
}

unless (defined $hmm_path) {
    ## Set a default path for HMMER. Since we want to build global models, set to hmmer-2.3.2.
    ## DNA/DNA models are available through 'nhmmer' in hmmer-3, but global models are not expected until hmmer-4 is released.
    $hmm_path = "/usr/local/hmmer/hmmer-2.3.2";
}

my $hmmer_version = find_hmmer($hmm_path);

if ($hmmer_version > 2.3 && defined $global) {                                            
    # Minor release numbers like 2.3.2 are not numeric, so we just look to see if 
    # one is trying to construct a global model with HMMER 3.0, which is not implemented.
    print "\nERROR: Cannot create global models with HMMER $hmmer_version. Exiting.\n\n"; 
    usage();
    exit(1);
}

my $date = POSIX::strftime("%m_%d_%Y_%H_%M_%S", localtime);

if ($clean && $header) {
    print "\nERROR: Options --header and --clean can not be used in conjuction. Exiting.\n\n";
    usage();
    exit(1);
}

my ($LTR_aln_files)  = get_files($in_dir);     ## use File::Spec to create new file paths 
my ($fas_files) = get_files($fas_dir);

if (scalar(@$LTR_aln_files) < 1 && scalar(@$fas_files) < 1) {
    print "\nERROR: No files were found in the input directories. Exiting.";
    usage();
    exit(1);
}

#chdir($in_dir); ## why not create a new dir for results?
my $aln_stats = get_aln_len(@LTR_aln_files); # return a hash-ref
    

# make one directory
my $model_dir = "LTR_hmmer2_models_".$date;
unless ( -e $model_dir ) {
    mkdir $model_dir || die "ERROR: Could not make directory: $model_dir\n";
}

chdir $model_dir;

for my $LTR_aln (@LTR_aln_files) {
    
    copy("../$LTR_aln",".") || die "Copy failed: $!";
    
    if (!defined $quiet) {
	print "\n\n============= building model for $LTR_aln ...\n";
    }
    build_model($LTR_aln,$hmm_path);
    
    unlink $LTR_aln;	
}

opendir my $hmms, "../$model_dir" || die "ERROR: Could not open directory: $model_dir\n";
my @LTR_hmm_files = grep /\.hmm.*$/, readdir $hmms;       # this is safe because models go into a specific directory based on 
closedir $hmms;                                           # execution time

my $hmmresults = "hmmsearch_out_".$date;
unless ( -e $hmmresults ) {
    mkdir $hmmresults;
}
chdir $hmmresults;

my $hmmresults_parsed = "hmmsearch_out_parsed_".$date;
unless ( -e $hmmresults_parsed ) {
    mkdir $hmmresults_parsed;
}

for my $hmm (@LTR_hmm_files) {
   
    my $indiv_results_dir = $hmm;
    $indiv_results_dir =~ s/\.hmm.*$/\_search_out/;
    mkdir($indiv_results_dir);
    chdir($indiv_results_dir);
    
    for my $fas (@fas_files) {
	
	my $hmmsearch_out = $hmm;
	$hmmsearch_out =~ s/\.hmm.*$//;
	$hmmsearch_out .= "-".$fas;
	$hmmsearch_out =~ s/\.fasta$/\.hmmer/;
	
	search_with_models($hmmsearch_out,$hmm,$hmm_path,$fas,$fas_dir);
	
	if ($report) {   
	    
	    if (!defined $quiet) {
		print "\n============= parsing report $hmmsearch_out ...\n";
	    }
	    
	    my $parsed = $hmmsearch_out;
	    $parsed =~ s/\.hmmer$/\_hmmer_parsed.txt/;
	    
	    my $seqfile = $parsed;
	    $seqfile =~ s/\.txt$/\_seq.fasta/;
	    
	    if ($hmm =~ m/\.hmms$/) {
		my $model_type = "global";
		write_hmmsearch_report($seqfile,$hmmsearch_out,$hmmresults_parsed,$model_type);
	    }
	    elsif ($hmm =~ m/\.hmmfs$/) {
		my $model_type = "multi-hit_local";
		write_hmmsearch_report($seqfile,$hmmsearch_out,$hmmresults_parsed,$model_type);
	    }
	    elsif ($hmm =~ m/\.hmmls$/) {
		my $model_type = "local";
		write_hmmsearch_report($seqfile,$hmmsearch_out,$hmmresults_parsed,$model_type);
	    }
	    
	    #write_hmmsearch_report($seqfile,$hmmsearch_out,$hmmresults_parsed,$model_type);
	    
	} next;
	
    }
    chdir "..";      # This is why not searching all seqs before
    
}

if ($clean) {
    my $parsed_empty_count = `find $hmmresults_parsed -size 0c | wc -l`;
    chomp $parsed_empty_count;
    
    if ($parsed_empty_count > 0) {
	system("find $hmmresults_parsed -size 0c | xargs rm");
    }
}

if ($report) {
    my $hmmsearch_summary = "all_hmmer_parsed_summary.out";
    #system("cat $hmmresults_parsed/*txt > $hmmresults_parsed/$hmmsearch_summary");
    chdir($hmmresults_parsed);
    system("ls | xargs cat > $hmmsearch_summary");
}

exit;

#
# methods
#
sub find_hmmer {
    my $hmm_path = shift;

    unless ($hmm_path =~ /\/$/ ) {
        $hmm_path = $hmm_path."/";
    }
    my $release;
    my $version = `$hmm_path/hmmbuild -h | head -2 | tail -1`;

    if ($version =~ /# HMMER (\d\.\d.*) \(/) {  # version 3.0
	$release = substr($1,0,3);                    
	return $release;
    } 
    elsif ($version =~ /HMMER (\d\.\d.*) \(/) { # version 2.3x
	$release = substr($1,0,3);
	return $release;
    } else {
	print "\nERROR: Could not determine HMMER Version. Check input path. Exiting.\n\n";
	usage();
	exit(1);
    }
}

sub get_files {
    my $indir = shift;
    opendir my $dir, $indir ); # || print "Can't open directory:\n$indir\n";
    my @dir_contents = grep /\.aln$|\.fa$|\.fas|\.fasta$/, readdir $dir;
    closedir $dir;

    return(\@dir_contents);
}

sub get_aln_len {
    my @aln_files = @_;
    
    #my $aln_report = "all_aln_stats_$date.txt";
    my %aln_stats;
    #open( REPORT, ">>$aln_report" ) || die "Could not open file: $aln_report\n";

    for my $aligned (@aln_files) {

	my $ltr_retro = $aligned;
	$ltr_retro =~ s/\.aln$//;
	my $aln_in = Bio::AlignIO->new(-file   => $aligned,
				       -format => 'clustalw');

	while ( my $aln = $aln_in->next_aln() ) {

	    #my @alnmts = join("\t",($aligned,$aln->length));
	
	    $aln_stats{$ltr_retro} = $aln->length;
	    # BAC17_sub1_c2_c1666_12133_25902_LTRseqs.aln    307
	    # BAC17_sub1_c2_c1666_31060_43611_LTRseqs.aln    577
	    # BAC17_sub1_c4_c1666_15251_26091_LTRseqs.aln    2160
	    # BAC17_sub1_c4_c1666_2529_8183_LTRseqs.aln      494
	    # foreach(@alnmts) {
	    #   print REPORT $_."\n";
	    # }

	}
	
    }
    #close(REPORT);
    return(\%aln_stats);
}
	
sub build_model {
    my ($aln, $hmmer_dir) = @_;
 
    my $hmmname      = $aln;
    my $hmmbuild     = $hmmer_dir."hmmbuild";
    my $hmmcalibrate = $hmmer_dir."hmmcalibrate";
    #my $calibrate_cmd = "$hmmcalibrate $hmmname 2>&1 > /dev/null";
    
    if (defined $global) {
	
	$hmmname =~ s/\.aln$/\.hmms/;
	my $hmm_cmd = "$hmmbuild -g $hmmname $aln 2>&1 > /dev/null";
	system($hmm_cmd);
	
	if (defined $calibrate) {
	    if (!defined $quiet) {
		print "\n============= calibrating model $hmmname ...\n";
	    }
	    my $calibrate_cmd = "$hmmcalibrate $hmmname 2>&1 > /dev/null";
	    system($calibrate_cmd);
	}
	
    } 
    elsif (defined $multi) {

	$hmmname =~ s/\.aln$/\.hmmfs/;
	my $hmm_cmd = "$hmmbuild -f $hmmname $aln 2>&1 > /dev/null";
	system($hmm_cmd);

	if (defined $calibrate) {
	    if (!defined $quiet) {
		print "\n============= calibrating model $hmmname ...\n";
	    }
	    my $calibrate_cmd = "$hmmcalibrate $hmmname 2>&1 > /dev/null";
	    system($calibrate_cmd);
	}
    } else {
	
	$hmmname =~ s/\.aln$/\.hmmls/;
	my $hmm_cmd = "$hmmbuild $hmmname $aln 2>&1 > /dev/null";
	system($hmm_cmd);
	
	if (defined $calibrate) {
	    if (!defined $quiet) {
		print "\n============= calibrating model $hmmname ...\n";
	    }
	    my $calibrate_cmd = "$hmmcalibrate $hmmname 2>&1 > /dev/null";
	    system($calibrate_cmd);
	}

    }
}

sub search_with_models {
    my ($search_out,$hmm_model,$hmmer_dir,$fasta,$fasta_dir) = @_;

    # hmmsearch
    if (!defined $quiet) {
	print "\n\n============= searching $fasta with model $hmm_model ...\n";
    }
    my $hmmsearch = $hmmer_dir."hmmsearch";

    my $hmmsearch_cmd = "$hmmsearch ../../$hmm_model ../../../../$fasta_dir$fasta > $search_out"; 
    system($hmmsearch_cmd);
    
    return $search_out;
}

sub write_hmmsearch_report {
    my ($seq_out,$search_report,$search_results_dir,$search_type) = @_;
    
    my $parsed = $search_report;
    $parsed =~ s/\.hmmer$/\_hmmer_parsed.txt/;

    my $ID_length_report = "all_aln_stats_$date.txt";

    my $match_threshold = defined($percentID_threshold) ? $percentID_threshold : "0.39";
   
    open my $out, ">$parsed" or die "\nERROR: Could not open file: $parsed\n";
    
    if (defined $seq) {
	print "\n============= writing sequence $seq_out ...\n" unless $quiet;
	open my $seq, ">$seq_out" or die "\nERROR: Could not open file: $seq_out";
    }

    my $hmmerin = Bio::SearchIO->new(-file => "$search_report", -format => 'hmmer');

    if (defined $header) {
	say $out join "\t", "#query", "query_length", "number_of_hits", "hit_name", 
	"perc_coverage","hsp_length", "hsp_perc_ident","hsp_query_start", "hsp_query_end", 
	"hsp_hit_start", "hsp_hit_end","search_type";
    }
    
    open my $alnin, "../../../$ID_length_report" or die "\nERROR: Could not open file: $ID_length_report\n";

    my $positions = 0;
	
    while( my $result = $hmmerin->next_result() ) {
	
	my $query      = $result->query_name();
	my $num_hits   = $result->num_hits();
	
	while( my $hit = $result->next_hit() ) {
	    
	    my $hitid    = $hit->name();
	    #my $signif   = $hit->significance();
	    
	    while( my $hsp = $hit->next_hsp() ) {
		
		#my $hsppid    = $hsp->percent_identity();
		#my $hspeval   = $hsp->evalue();                   # appears to be redundant with $signif
		#my $hspident  = $hsp->frac_identical;
		my @ident_pos = $hsp->seq_inds('query','identical');
		
		my $hspgaps   = $hsp->gaps;
		my $hsplen    = $hsp->length('total');
		my $hstart    = $hsp->start('hit');
		my $hstop     = $hsp->end('hit');
		my $qstart    = $hsp->start('query');
		my $qstop     = $hsp->end('query');
		my $qstring   = $hsp->query_string;
		
		while (<$alnin>) {
		    chomp;
		    my @id_len_from_aln = split("\t",$_);
		    my $aln_id = $id_len_from_aln[0];
		    $aln_id =~ s/\.aln$//;
		    my $aln_len = $id_len_from_aln[1];
		    my $percent_coverage = sprintf("%.2f",$hsplen/$aln_len);
		    
		    if ($aln_id =~ m/$query/) {
			
			$positions++ for @ident_pos;
			
			if ( ($hsplen >= 80) && ($hsplen >= $aln_len * 0.80)  && ($hsplen == 0) ) { ## Need to revisit this $hsplen == 0 business
			    
			    my $percent_identity = 0;
			    
			    say $out join "\t", $query, $aln_len, $num_hits, $hitid, $percent_coverage,
			    $hsplen, $percent_identity,$qstart, $qstop, $hstart, $hstop,$search_type;
			    
			} else {
			    if ( ($hsplen >= 80) && ($hsplen >= $aln_len * 0.80) ) {
				#my $percent_identity = sprintf("%.2f",($hsplen - $hspgaps)/$hsplen);
				my $percent_identity = sprintf("%.2f",$positions/$hsplen);
				
				if ($percent_identity >= $match_threshold) {
				    
				    say $out join "\t", $query, $aln_len, $num_hits, $hitid, $percent_coverage,
				    $hsplen, $percent_identity,$qstart, $qstop, $hstart, $hstop,$search_type;

				    if (defined $seq) {
					my $seqid = ">".$query."|".$hitid."_".$hstart."-".$hstop; ## It makes more sense to show the location of the hit
                                                                                                  ## Also, this would pave the way for creating a gff of solo-LTRs
                                                                                                  ## my $seqid = ">".$query."|".$hitid."_".$hstart."-".$hstop
 					say $seq join "\n", $seqid, $qstring;
				    }

				    
				}
				
				# don't want to print every sequence, only the matches!!!!!!!
                                #if (defined $seq) {
				#    my $seqid = ">".$query."|".$hitid."_".$qstart."-".$qstop; ## It makes more sense to show the location of the hit 
                                                                                              ## Also, this would pave the way for creating a gff of solo-LTRs
                                                                                              ## my $seqid = ">".$query."|".$hitid."_".$hstart."-".$hstop
				    
				#    print SEQ join("\n", ($seqid, $qstring)),"\n";
				    
				#}
			    }
			}
		    }
		}
	    }
	}
    }
    close $alnin;
    close $out;
    close $seq;
    
    if ( -s $parsed ) {
	copy("$parsed","../$search_results_dir") || die "Copy failed: $!";
    } else {
	unlink $parsed if defined $clean;
    }
    if (defined $seq) {	
	my $seq_dir = "hmmer_domain_sequences-".$date;
	my @ext_domain_seqs = glob "*fasta";

	if ( -d "../$seq_dir" ) {	   	    
	    for my $ext_seq (@ext_domain_seqs) {
		if ( -s $ext_seq ) {
		    copy("$ext_seq","../$seq_dir") || die "Copy failed: $!"; 
		} else {
		    if (defined $clean) {
			unlink($ext_seq);
		    }
		}
	    }
	} else {
	    system("mkdir ../$seq_dir");
	    
	     for my $ext_seq (@ext_domain_seqs) {
		 if ( -s $ext_seq ) {
		     copy("$ext_seq","../$seq_dir") || die "Copy failed: $!";
		 } else {
		     unlink $ext_seq if defined $clean;
		 }
	     }
	}
	if (defined $clean) {
	    my $domain_seqs_empty_count = `find ../$seq_dir -size 0c | wc -l`;
	    chomp $domain_seqs_empty_count;

	    if ($domain_seqs_empty_count > 0) {		
		system("find ../$seq_dir -size 0c | xargs rm");
	    }
	}
    }
}

sub usage {
    my $script = basename($0);
    print STDERR <<END
USAGE: $script -i /path/to/dir/of/LTR/alns -f /path/to/dir/of/fastas/to/search <options>

Required:
   -i|indir         :    The input directory of LTR alignments in clustalw format.
                         Must have *aln ending.
   -f|fasdir        :    Input directory of fasta files to search.
                         Must have *fas or *fasta ending.

Options:
   -hmm_path             Path to the directory containing the HMMER binaries			 
   -p|percent_ident :    Percent identity threshold for matches. (default 0.39).
                         NB: For a threshold of 80 percent say 0.80.
   --quiet          :    Be quiet, don't print all of the progress to screen.
   --report         :    Parse hmmsearch of each sequence and produce a summary of align statistics.
   --global         :    Search with model of global alignments w.r.t the sequence and hmm. Will match
                         complete domains.
   --multi          :    Do multi-hit local alignment. Good for identifying fragments.
   --seq            :    Extract query sequence from domain alignment.
   --header         :    Print a column header in parsed hmmsearch report. Cannot be used with --clean option.
   --clean          :    Remove 0 byte files. This limits all of the empty files created from negative search
                         results. Will not work with the --header option though.
   --calibrate      :    Calibrate each model prior to searching.

END
}
