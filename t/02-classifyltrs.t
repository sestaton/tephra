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

use Test::More tests => 3;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd      = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir  = File::Spec->catdir('t', 'test_data');
my $outdir   = File::Spec->catfile($testdir, 't_family_domains');
my $genome   = File::Spec->catfile($testdir, 'ref.fas');
my $ingff    = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered.gff3');
my $outgff   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified.gff3');
my $outfas   = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified.fasta');
my $repeatdb = File::Spec->catfile($testdir, 'repdb.fas');

SKIP: {
    skip 'skip development tests', 3 unless $devtests;
    my @assemb_results = capture { system([0..5], "$cmd classifyltrs -h") };

    ok(@assemb_results, 'Can execute classifyltrs subcommand');

    my $find_cmd = "$cmd classifyltrs -g $genome -d $repeatdb -i $ingff -o $outgff -r $outdir";
    #say STDERR $find_cmd;
    
    my @ret = capture { system([0..5], $find_cmd) };

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
    
    ok( $ct == 6, 'Correct number of classified elements in combined family file' );
    say "$ct total combined elements";

    unlink $outfas;

};

unlink $ingff;
done_testing();
