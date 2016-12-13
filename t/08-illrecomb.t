#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Path          qw(remove_tree);
use File::Find;
use File::Spec;
use File::Copy;
#use Data::Dump::Color;

use Test::More tests => 8;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $outdir    = File::Spec->catdir($testdir,  't_family_domains');
my $allstfile = File::Spec->catfile($testdir, 'gypsy_illrecomb_stats.tsv');
my $illstfile = File::Spec->catfile($testdir, 'gypsy_illrecomb_illrecstats.tsv');
my $seqfile   = File::Spec->catfile($testdir, 'gypsy_illrecomb_seqs.fasta');
my $genome    = File::Spec->catfile($testdir, 'ref.fas');
my $testfile  = File::Spec->catfile($testdir, 'tephra_ltrs_gypsy_family9.fasta');

SKIP: {
    skip 'skip development tests', 8 unless $devtests;

    my @results = capture { system([0..5], "$cmd illrecomb -h") };    
    ok(@results, 'Can execute illrecomb subcommand');

    my $find_cmd = "$cmd illrecomb -i $testfile -s $allstfile -r $illstfile -o $seqfile";
    #say STDERR $find_cmd;

    my @ret = capture { system([0..5], $find_cmd) };
    #system([0..5], $find_cmd);

    ok( -s $allstfile, 'Generated statistics for all gap sites' );
    ok( -s $illstfile, 'Generated statistics for all putative illegetimate recombination sites' );
    ok( -s $seqfile,   'Generated sequences flanking all putative illegetimate recombination sites' );
    
    my $seqct = 0;
    open my $in, '<', $seqfile;
    while (<$in>) { $seqct++ if /^>/; }
    #say STDERR "seqct: $seqct";
    ok( $seqct == 36, 'Correct number of illigetimate recombination events detected' );
    close $in;
    
    my ($qmatch, $hmatch) = (0, 0);
    open my $stats, '<', $illstfile;
    while (<$stats>) {
	chomp;
	$qmatch++ if /^Query match string/;
	$hmatch++ if /^Hit match string/;
    }

    #say STDERR join q{ }, $qmatch, $hmatch;
    ok( $qmatch == 18, 'Correct number of illigetimate recombination events detected upstream of gap' );
    ok( $hmatch == 18, 'Correct number of illigetimate recombination events detected downstream of gap' );
    ok( $seqct == $qmatch+$hmatch, 'Correct number of illigetimate recombination events detected' );
    
    unlink $allstfile, $illstfile, $seqfile;

    ## clean up
    my @outfiles;
    find( sub { 
        push @outfiles, $File::Find::name 
            if /\.log$|RLG_family9/ }, $testdir);

    unlink @outfiles;
};

done_testing();
