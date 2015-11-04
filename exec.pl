#!/home/statonse/perl5/perlbrew/perls/perl-5.22.0/bin/perl

use strict; 
use warnings;
use Config;
use Bio::SearchIO;
use File::Spec;
use File::Basename;
use List::Util          qw(min);
use File::Copy          qw(copy);
use IPC::System::Simple qw(system capture);
use Try::Tiny;
use Getopt::Long;
use Carp 'croak';

my $seq = 'TGAATGTGCTCTCAGTTAACATTCGAGGTTTGGGTCGGGTTGATAAGGGGGATTGGATTTCTAATATTCGGGTGAAAAACGAGGCGTCTTTCGTTATGTTACAAGAGACTCAATTTGCTTCGTTGCAAGGTGTGGATATTGGTAAATATTGGGGTAGCGGTACGTTTGATTCTGAGTTTATGAATGCCATGGGTCGGTCCGGGGGTCTTCTATCTCTTTGGGATCAAAAGTTGTTTCAGAAGAGTTCGGTGCAGAAACATAGATACTATCTCGCTGTGCATGGGTATATTAAAGGTAACGGTGCGAAGGTGTGTTTGGTGAATGTATACGCTCCTCAAAAAATTACAGAAAAAAGACTATTGTGGCGGGAATTGGAAAGGCTTGTGCAGCAAGACGAGACATATTGGGTAGTGGGCGGTGATTTTAATTGTGTCCGGGATAGAAGTGAGAGGAAGAATACAAACTTCAACGCCGCGTCTTCAATTGAATTTAATGATTTCTTAGATGGGGTTGGTTTGCATGAGTTCGGTTTA';

my $phmm_file = 'MGEScan_nonLTR_v2/my_rewrite2/pHMM/RTE.en.hmm';
my $pep_file  = 'testtr.faa';
my $seq_file  = 'testtr.fas';
open my $o, '>', $seq_file;
say $o $seq;
close $o;

my $exitv;
try {
    $exitv = system([0..5], "transeq", "-frame=f", $seq_file, "-outseq=$pep_file", "-auto");
}
catch {
    print "\nERROR: transeq failed with exit value: $exitv. Here is the exception: $_\n";
    exit;
};

my $hmmsearch   = _find_hmmsearch();
my @hmm_results = capture([0..5], $hmmsearch, $phmm_file, $pep_file);
_parse_hmmsearch(\@hmm_results);


sub _parse_hmmsearch {
    my ($hmm_results) = @_;

    my $hmmout = 'testout.txt';
    open my $o, '>', $hmmout;
    for my $res (@$hmm_results) {
	print $o $res;
    }
    close $o;

    my $hmmer_in = Bio::SearchIO->new(-file => $hmmout, -format => 'hmmer');

    my @evalues;
    while ( my $result = $hmmer_in->next_result ) {    
	while ( my $hit = $result->next_hit ) {
	    my $score  = $hit->raw_score;
	    my $signif = $hit->significance;
	    while ( my $hsp = $hit->next_hsp ) {
		my $e_val = $hsp->evalue;
		push @evalues, $e_val;
	    }
	}
    }
    my $best_hit = min(@evalues);
    print $best_hit;
}

sub _find_hmmsearch {
    my $hmmer2_env = $ENV{HMMER2} // '/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2';
    my $hmmsearch = join "/", $hmmer2_env, 'src', 'hmmsearch';
    if (-e $hmmsearch && -x $hmmsearch) {
	return $hmmsearch;
    }
    else {
	$hmmsearch = join "/", $hmmer2_env, 'bin', 'hmmsearch';
	if (-e $hmmsearch && -x $hmmsearch) {
	    return $hmmsearch;
	}
	else {
	    croak "\nERROR: Could not find 'hmmsearch'. Make sure variable 'HMMER2' points to the hmmer2 root".
		" directory. Exiting.\n";
	}
    }
}
