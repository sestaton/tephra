#!/home/statonse/perl5/perlbrew/perls/perl-5.22.0/bin/perl

use 5.010;
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

my $hmmsearch   = _find_hmmsearch();

#my @hmm_results = capture([0..5], $hmmsearch, $phmm_file, $pep_file);
my $file = shift;
open my $in, '<', $file;
my @hmm_results = <$in>;
close $in;
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
	my $query = $result->query_name;
	while ( my $hit = $result->next_hit ) {
	    my $hitid  = $hit->name;
	    my $score  = $hit->raw_score;
	    my $signif = $hit->significance;
	    my $bits   = $hit->bits;
	    while ( my $hsp = $hit->next_hsp ) {
		my $hstart  = $hsp->start('hit');
		my $hstop   = $hsp->end('hit');
		my $qstart  = $hsp->start('query');
		my $qstop   = $hsp->end('query');
		my $score   = $hsp->score;
		my $e_val = $hsp->evalue;
		#push @evalues, $e_val;
		say join q{ }, $hitid, $hstart, $hstop, $qstart, $qstop, $score, $e_val;
	    }
	}
    }
    #my $best_hit = min(@evalues);
    #print $best_hit;
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
