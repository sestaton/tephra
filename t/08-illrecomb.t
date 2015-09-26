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
use Data::Dump;

use Test::More tests => 1;

my $bindir = File::Spec->catdir('t', 'gt', 'bin');
local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my @results = capture { system([0..5], "$cmd illrecomb -h") };

ok(@results, 'Can execute illrecomb subcommand');

done_testing();
