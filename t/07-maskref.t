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

use Test::More tests => 4;

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
#my $repeatdb = File::Spec->catfile($testdir, 'repdb.fas');
my $repeatdb = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified.fasta');
my $masked   = File::Spec->catfile($testdir, 'ref_masked.fas');
my $log      = File::Spec->catfile($testdir, 'ref_masked.fas.log');

{
    my @help_args = ($cmd, 'maskref', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute maskref subcommand');
}

my @mask_cmd = ($cmd, 'maskref', '-g', $genome, '-d', $repeatdb, '-o', $masked, $log);
say STDERR join q{ }, @mask_cmd if $devtests;
my ($stdout, $stderr, $exit) = capture { system([0..5], @mask_cmd) };

ok( -e $masked, 'Can mask reference' );
ok( -e $log,    'Can mask reference' );

for my $line (split /^/, $stdout) { 
    if ($line =~ /Total masked bases:  (\d+.\d+)%/) {
	ok( $1 > 20, 'Logged non-zero masked bases of output masked genome' );
    }
}

## clean up
unlink $log, $repeatdb;

done_testing();
