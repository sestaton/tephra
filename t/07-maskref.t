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

use Test::More tests => 3;

$| = 1;

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my $repeatdb = File::Spec->catfile($testdir, 'repdb.fas');
my $masked   = File::Spec->catfile($testdir, 'ref_masked.fas');
my $log      = File::Spec->catfile($testdir, 'ref_masked.fas.log');

{
    my @help_args = ($cmd, 'maskref', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute maskref subcommand');
}

my @mask_cmd = ($cmd, 'maskref', '-g', $genome, '-d', $repeatdb, '-o', $masked, $log);
#say STDERR $mask_cmd;
my @ret = capture { system([0..5], @mask_cmd) };

ok( -e $masked, 'Can mask reference' );
ok( -e $log, 'Can mask reference' );

## clean up
unlink $log;

done_testing();
