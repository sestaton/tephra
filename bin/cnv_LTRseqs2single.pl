#!/usr/bin/env perl

## combine with cnv_LTRseqs2align.pl script

use strict;
use warnings;
use Bio::SeqIO;
#use Getopt::Long;  # get options such as: --align, --combine, --report
                    # default would be: -l -r or something similar


my $usage = "\nUSAGE: LTRseqs2align2.pl 5primeseqs.fasta 3primeseqs.fasta

\tThe files must be formatted with seq_add_name.pl to work with clustalw2.\n";

my $left = shift or die $usage;
my $right = shift or die $usage;


my %seq_in = ( 
               '5prime' => Bio::SeqIO->new('-file'   => "<$left",
                                           '-format' => 'fasta'),
               '3prime' => Bio::SeqIO->new('-file'   => "<$right",
					   '-format' => 'fasta'),
             );

while (my $fiveltr = $seq_in{'5prime'}->next_seq) {
     my $fiveid = $fiveltr->id;
     $fiveid =~ s/5prime_//;
     my $LTRseqs_name = $fiveid;
     $LTRseqs_name .= "_LTRseqs.fasta";
     my $seq_out = Bio::SeqIO->new(-file => ">$LTRseqs_name", -format => 'fasta');    

     while (my $threeltr = $seq_in{'3prime'}->next_seq) {
	 my $threeid = $threeltr->id;
	 $threeid =~ s/3prime_//;
	 if ($fiveid =~ /$threeid/) {		 
	     $seq_out->write_seq($fiveltr,$threeltr);
	     last;
	 }
     }
}
