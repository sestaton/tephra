#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Find;
use File::Spec;

use Test::More tests => 2;

my $bindir = File::Spec->catdir('t', 'gt', 'bin');
local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $model   = File::Spec->catfile($testdir, 'te.hmm');
my $trnas   = File::Spec->catfile($testdir, 'trnas.fas');

my @assemb_results = capture { system([0..5], "$cmd findtrims -h") };

ok(@assemb_results, 'Can execute findtrims subcommand');

my $find_cmd = "$cmd findtrims -g $genome -t $trnas -d $model --clean";
say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };

my @files;
find( sub { push @files, $File::Find::name if /trim/ }, $testdir);
ok( @files == 6, 'Can find some trims' ); # 2 ltrdigest files + combined file

## clean up
my @outfiles;
find( sub { push @outfiles, $File::Find::name if /^ref_trim/ }, $testdir);
unlink @outfiles;
    
done_testing();
