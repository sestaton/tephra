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
#use Data::Dump::Color;

use Test::More tests => 5;

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my $repeatdb = File::Spec->catfile($testdir, 'repdb.fas');
my $masked   = File::Spec->catfile($testdir, 'ref_masked.fas');
my $outfile  = File::Spec->catfile($testdir, 'ref_masked_tephra_repdb_fragments.gff3');
my $log      = File::Spec->catfile($testdir, 'ref_findfragments.log');
my $thrlog   = File::Spec->catfile($testdir, 'tephra_fragment_searches.log');

SKIP: {
    skip 'skip development tests', 5 unless $devtests;
    {
        my @help_args = ($cmd, 'findfragments', '-h');
        my ($stdout, $stderr, $exit) = capture { system(@help_args) };
        #say STDERR "stderr: $stderr";
        ok($stderr, 'Can execute findfragments subcommand');
    }

    my $frag_cmd = "$cmd findfragments -g $masked -d $repeatdb -o $outfile > $log";
    ##say STDERR $frag_cmd;    
    my @ret = capture { system([0..5], $frag_cmd) };

    ok( -e $outfile, 'Can run findfragments and generate GFF3 output' );
    ok( -e $log,     'Expected log for findfragments command produced' );
    ok( -e $thrlog,  'Expected log for findfragments command produced' );
    
    open my $in, '<', $outfile;
    my $fragct = 0;
    while (<$in>) {
	chomp;
	next if /^#/;
	$fragct++;
    }

    ok( $fragct == 3, 'Expected number of fragments discovered' );

    ## clean up
    unlink $log, $thrlog, $masked, $outfile;
};

done_testing();
