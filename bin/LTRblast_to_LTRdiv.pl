#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Getopt::Long;

my $infile;
my $outfile;
my $usage = "USAGE: $0 -i inreport -o parsedreport

\tFirst, run the script seq_add_name.pl on the file of 5prime sequences
\tand 3prime sequences separately, including the name 5prime and 3prime, respectively. 
\tThis will make it possible to distinguish the five-prime and three-prime LTR sequences 
\tin the blast report (and ignore other hits). 

\tNext, BLAST the LTR sequences together with a low e-value (e.g. 1e-20) and include
\tthe -m 8 option when running blastall. 

\tLast, simply run this script on the BLAST report.";

GetOptions(
           'i|infile=s'   => \$infile,
           'o|outfile=s'  => \$outfile,
          );

# open the infile or die with a usage statement
if (!$infile || !$outfile) {
    say "\nERROR: Command line not parsed correctly. Exiting.\n";
    say $usage;
    exit(1);
}

open my $in, '<', $infile or die "\nERROR: Could not open file: $infile.\n";
open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile.\n";

say $out join "\t", "#5primeLTR", "3primeLTR", "Percent_identity", "Alignment_length";

while (<$in>) { 
    chomp; 
    my @blfields = split /\t/;

    my $fiveprime = $blfields[0];
    my $threeprime = $blfields[1];

    if ($blfields[0] !~ $blfields[1]) {
	$fiveprime =~ s/5prime$//;
	$threeprime =~ s/3prime$//;

        if ($fiveprime eq $threeprime) { 
	    say $out @blfields[0..3];
	 }
    }
}

exit;
