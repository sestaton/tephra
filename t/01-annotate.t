#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(capture);
use File::Path          qw(remove_tree);
use File::Find;
use File::Spec;

use Test::More tests => 14;

my $tot       = 0;
my $scr       = 0;
my $seqnum    = 20;
my $outdir    = "VelvetOpt_k59-k59";
my $pairfile  = File::Spec->catfile('t', 'test_data', 't_cpseqs_screened_paired_interl.fas');
my $upairfile = File::Spec->catfile('t', 'test_data', 't_cpseqs_screened_unpaired.fas');

my $cmd = File::Spec->catfile('bin', 'chloro');
my @assemb_results = capture([0..5], "$cmd assemble -p $pairfile -u $upairfile -s 59 -e 59");

my @log;
find( sub {
    push @log, $File::Find::name if -f and /logfile.txt$/;
      }, $outdir);

my ($hashv, $mapsize, $contignum, $n50, $longest, $total, $over1knum, $over1kbases);

for my $logfile (@log) {
    open my $in, '<', $logfile;
    
    while (<$in>) {
	chomp;
	if (/Velvet hash value: (\d+)/) {
	    $hashv = $1;
	}
	if (/Roadmap file size: (\d+)/) {
	    $mapsize = $1;
	}
	if (/Total number of contigs: (\d+)/) {
	    $contignum = $1;
	}
	if (/n50: (\d+)/) {
	    $n50 = $1;
	}
	if (/length of longest contig: (\d+)/) {
	    $longest = $1;
	}
	if (/Total bases in contigs: (\d+)/) {
	    $total = $1;
	}
	if (/Number of contigs > 1k: (\d+)/) {
	    $over1knum = $1;
	}
	if (/Total bases in contigs > 1k: (\d+)/) {
	    $over1kbases = $1;
	}
    }
    close $in;
}

is( $hashv,       59,     'Correct hash value used in assembly' );
is( $mapsize,     299868, 'Expected Roadmap file generated' );
is( $contignum,   128,    'Expected number of contigs produced' );
is( $n50,         1337,   'Expected N50 contig size produced' );
is( $longest,     5054,   'Expected size of longest contig' );
is( $total,       101938, 'Expected total bases in contigs' );
is( $over1knum,   31,     'Expected number of contigs >1kb' );
is( $over1kbases, 60362,  'Expected number of bases in contigs over 1kb');

my $vellog   = File::Spec->catfile($outdir, 'Log');
my $pregraph = File::Spec->catfile($outdir, 'PreGraph');
my $graph    = File::Spec->catfile($outdir, 'Graph');
my $graph2   = File::Spec->catfile($outdir, 'Graph2');
my $stats    = File::Spec->catfile($outdir, 'stats.txt');
my $contigs  = File::Spec->catfile($outdir, 'contigs.fa');

ok( -s $vellog,   'Created velvet log' );
ok( -s $pregraph, 'Created velvet pregraph' );
ok( -s $graph,    'Created velvet graph' );
ok( -s $graph2,   'Created velvet graph2' );
ok( -s $stats,    'Created velvet stats' );
ok( -s $contigs,  'Created velvet contigs from assembly' );

unlink $pairfile;
unlink $upairfile;
remove_tree($outdir);

done_testing();
