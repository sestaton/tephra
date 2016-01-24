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

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $model   = File::Spec->catfile($testdir, 'te.hmm');

my @results = capture { system([0..5], "$cmd findtirs -h") };

ok(@results, 'Can execute findtirs subcommand');

my $find_cmd = "$cmd findtirs -g $genome -d $model --clean";
#say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };

my @files;
find( sub { push @files, $File::Find::name if /tirs_?(?:filtered)?.gff3$/ }, $testdir);
ok( @files == 2, 'Can find some tirs' ); # original + filtered gff3

## clean up
my @outfiles;
find( sub { push @outfiles, $File::Find::name if /^ref_tirs/ && ! /filtered.gff3$/ }, $testdir);
unlink @outfiles;
    
done_testing();
