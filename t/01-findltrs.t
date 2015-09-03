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

use Test::More tests => 2;

my $bindir = File::Spec->catdir('t', 'gt', 'bin');
local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd    = File::Spec->catfile('bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome = File::Spec->catfile('t', 'test_data', 'ref.fas');
my $model  = File::Spec->catfile('t', 'test_data', 'te.hmm');
my $trnas  = File::Spec->catfile('t', 'test_data', 'trnas.fas');

my @assemb_results = capture { system([0..5], "$cmd findltrs -h") };

ok(@assemb_results, 'Can execute findltrs subcommand');

my $find_cmd = "$cmd findltrs -g $genome -t $trnas -p $model -o filtered.gff3 --clean";
#say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };
#dd \@ret;

my @files;
find( sub { push @files, $File::Find::name if /\.gff3/ }, $testdir);
ok( @files == 4, 'Can find some ltrs' );

## clean up
my @outfiles;
find( sub { push @outfiles, $File::Find::name if /^ref_ltr/ }, $testdir);
unlink @outfiles;
unlink glob("*.gff3"); # TODO
    
done_testing();
