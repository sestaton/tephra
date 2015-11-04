#!/usr/bin/env perl

#TODO: Take multiple hmm names at the command line.
#      Make use of the (undocumented) Easel library shipped with HMMER,
#      which provides access to many parsing/extraction utilities.
use 5.010;
use strict;
use warnings;
use Bio::SearchIO;                                                  
use Getopt::Long;                                                   
use File::Basename;

my $infile; 
my $outfile;
my $seqfile;
my $hmmname;
my $hmmseq;
my $seq;

GetOptions(
	   "i|infile=s"      => \$infile,
	   "o|outfile=s"     => \$outfile,
	   "s|seq=s"         => \$seqfile,
	   "h|hmmname=s"     => \$hmmname,
	   );

# open the infile or die with a usage statement
usage() and exit(1) if !$infile;

#open my $out, ">", $outfile or die "\nERROR: Could not open file: $!\n";
#if ($seqfile) { open $seq, '>', $seqfile or die "\nERROR: Could not open file: $!\n"; }
#if ($hmmname) { open $hmmseq, '>', $hmmname or die "\nERROR: Could not open file: $!\n"; }

my $hmmer_in = Bio::SearchIO->new(-file => $infile, -format => 'hmmer');

#say $out join "\t", "query", "query_length", "number_of_hits", "hit_name", 
#    "hit_acc", "hit_score", "hit_significance", "hsp_length", "hsp_query_start", 
#    "hsp_query_end", "hsp_hit_start", "hsp_hit_end";

while ( my $result = $hmmer_in->next_result ) {    
    my $query    = $result->query_name;
    my $qlen     = $result->query_length;
    my $num_hits = $result->num_hits;
       
    while ( my $hit = $result->next_hit ) {
	my $hitid  = $hit->name;
	my $hitacc = $hit->accession;
	my $score  = $hit->score;
	my $signif = $hit->significance;
	
	while ( my $hsp = $hit->next_hsp ) {
	    my $hsplen  = $hsp->length('total');
	    my $hstart  = $hsp->start('hit');
	    my $hstop   = $hsp->end('hit');
	    my $qstart  = $hsp->start('query');
	    my $qstop   = $hsp->end('query');
	    my $qstring = $hsp->query_string;
	    my $hstring = $hsp->hit_string;
	    my $evalue  = $hsp->evalue;
	    #my $seqid = ">".$query."|".$hitid."_".$qstart."-".$qstop;
	     
	    say join "\t", $hitid, $hstart, $hstop;
		
	    #if ($seqfile) {
		#say $seq join "\n", $seqid, $qstring;
		
	    #}
	    #if ($hmmname) {
		#say $hmmseq join "\n", $seqid, $qstring
		#    if $hmmname =~ /$hitid/;
	    #}
	}
    }
}
#close $out;
#close $seq if $seqfile;
#close $hmmseq if $hmmname;

exit;
#
# methods
#
sub usage {
    my $script = basename($0);
    print STDERR <<END

USAGE: $script -i myresults.hmmer -o myparsedresults.txt [-s] [-h]

Required:
    -i|infile        :    HMMER report to parse.
    -o|outfile       :    File name to write the parsed results to.
    
Options:
    -s|seq           :    Write the matching query string to a separate Fasta file.
    -h|hmmname       :    A name of domain to search for in the results. If given,
                          all matches to the domain will be written to a file.

END
}

