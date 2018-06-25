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

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $outdir    = File::Spec->catdir($testdir,  'ltr_family_domains');
my $allstfile = File::Spec->catfile($testdir, 'gypsy_illrecomb_stats.tsv');
my $illstfile = File::Spec->catfile($testdir, 'gypsy_illrecomb_illrecstats.tsv');
my $seqfile   = File::Spec->catfile($testdir, 'gypsy_illrecomb_seqs.fasta');
my $genome    = File::Spec->catfile($testdir, 'ref.fas');
my $testfile  = File::Spec->catfile($testdir, 'tephra_ltrs_gypsy_family9.fasta');
my $log       = File::Spec->catfile($testdir, 'all_illrecomb_muscle_reports.log');

SKIP: {
    skip 'skip development tests', 8 unless $devtests;
    {
        my @help_args = ($cmd, 'illrecomb', '-h');
        my ($stdout, $stderr, $exit) = capture { system(@help_args) };
        #say STDERR "stderr: $stderr";
        ok($stderr, 'Can execute illrec subcommand');
    }

    my @find_cmd = ($cmd, 'illrecomb', '-i', $testfile, '-s', $allstfile, '-r', $illstfile, '-o', $seqfile);
    say STDERR join q{ }, @find_cmd;

    my @ret = capture { system([0..5], @find_cmd) };
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
    
    unlink $allstfile, $illstfile, $seqfile, $log;

    ## clean up (this is done automatically now in v0.06.0)
    #my @outfiles;
    #find( sub { 
    #    push @outfiles, $File::Find::name 
    #        if /\.log$|RLG_family9/ }, $testdir);
    
    #unlink @outfiles;
};

done_testing();
