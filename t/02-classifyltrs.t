#!/usr/bin/env perl

## todo: add test for exemplars

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

use Test::More tests => 4;

my $bindir = File::Spec->catdir('t', 'gt', 'bin');
local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $outdir   = File::Spec->catfile($testdir, 't_family_domains');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my $gff      = File::Spec->catfile($testdir, 'ref_ltrdigest85_combined_filtered.gff3');
my $repeatdb = File::Spec->catfile($testdir, 'repdb.fas');

my @assemb_results = capture { system([0..5], "$cmd classifyltrs -h") };

ok(@assemb_results, 'Can execute classifyltrs subcommand');

my $find_cmd = "$cmd classifyltrs -g $genome -d $repeatdb -f $gff -o $outdir";
#say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };

my @files;
find( sub { 
    push @files, $File::Find::name if /(?:gypsy|copia|unclassified).gff3$/ }, 
      $testdir );
ok( @files == 3, 'Correctly generated GFF3 files for LTR superfamilies' ); # 2 ltrdigest files + combined file

my @dirs;
find( sub { 
    push @dirs, $File::Find::name if -d && /(?:gypsy|copia|unclassified)$/ }, 
      $outdir );
ok( @files == 3, 'Correctly generated subdirectories for family level classification' );

my @all;
find( sub {
    push @all, $File::Find::name if /ref_combined_LTR_families.fasta/ },
      $outdir );

my $combined = shift @all;
open my $in, '<', $combined;
my $ct = 0;
while (<$in>) { $ct++ if /^>/; }
close $in;

ok( $ct == 6, 'Correct number of classified elements in combined family file' );

## clean up
my @outfiles;
find( sub { 
    push @outfiles, $File::Find::name 
	if /^ref_ltr/ && ! /(?:gypsy|copia|unclassified).gff3$/ && ! /$gff/ }, $testdir);

#unlink @outfiles;
#remove_tree( $outdir, { safe => 1 } );
    
done_testing();
