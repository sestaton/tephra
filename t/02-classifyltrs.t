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
#use Data::Dump::Color;

use Test::More tests => 9;

$| = 1;

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $outdir   = File::Spec->catfile($testdir, 't_family_domains');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my $log      = File::Spec->catfile($testdir, 'ref_tephra_classifyltrs.log');
my $ingff    = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered.gff3');
my $outgff   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified.gff3');
my $outfas   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified.fasta');
my $uncfas   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_unclassified.fasta');
my $copdom   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_copia_domain_org.tsv');
my $gypdom   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_gypsy_domain_org.tsv');
my $uncdom   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_unclassified_domain_org.tsv');
my $famdom   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified_family-level_domain_org.tsv');
my $repeatdb = File::Spec->catfile($testdir, 'repdb.fas');

{
    my @help_args = ($cmd, 'classifyltrs', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute classifyltrs subcommand');
}

my @find_cmd = ($cmd, 'classifyltrs', '-g', $genome, '-d', $repeatdb, '-i', $ingff, '-o', $outgff, '-r', $outdir);
#say STDERR join q{ }, @find_cmd;

my @ret = capture { system([0..5], @find_cmd) };

my @dirs;
find( sub { 
    push @dirs, $File::Find::name if -d && /(?:gypsy|copia|unclassified)$/ }, 
      $outdir );
ok( @dirs == 3, 'Correctly generated subdirectories for family level classification' );
#say scalar(@dirs)," number of subdirectories";


open my $in, '<', $outfas;
my $ct = 0;
while (<$in>) { $ct++ if /^>/; }
close $in;

open my $din, '<', $famdom;
my $fct = 0;
while (<$din>) { $fct++ if /\S+/; }
close $din;

#say STDERR "ct: $ct";
ok( $ct == 6, 'Correct number of classified elements in combined family file' );
ok( -e $log, 'Generated log for classifyltrs command' );
ok( -e $copdom, 'Generated domain organization file for Copia elements' );
ok( -e $gypdom, 'Generated domain organization file for Gypsy elements' );
ok( -e $uncdom, 'Generated domain organization file for Unclassfied elements' );
ok( -e $famdom, 'Generated domain organization file individual elements' );
ok( $ct == $fct, 'Correct number of classified elements in family-level domain organization file' );

unlink $outfas, $ingff, $log, $uncfas, $copdom, $gypdom, $uncdom, $famdom;

done_testing();
