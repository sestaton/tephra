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

use Test::More tests => 5;

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $outdir    = File::Spec->catdir($testdir, 't_family_domains');
my $resdir    = File::Spec->catdir($outdir, 'ref_ltrdigest85_combined_filtered_gypsy');
my $modeldir  = File::Spec->catdir($resdir, 'Tephra_LTR_exemplar_models');
my $allstfile = File::Spec->catfile($resdir, 'gypsy_sololtr_stats.tsv');
my $seqfile   = File::Spec->catfile($modeldir, 
				    'RLG_family0_exemplar_ltrs_clustal-out_ref_masked_hmmer_parsed_seq.fasta');
my $parsfile  = File::Spec->catfile($modeldir,
				    'RLG_family0_exemplar_ltrs_clustal-out_ref_masked_hmmer_parsed.txt');
my $masked    = File::Spec->catfile($testdir, 'ref_masked.fas');
my @results   = capture { system([0..5], "$cmd sololtr -h") };

ok(@results, 'Can execute sololtr subcommand');

my $find_cmd = "$cmd sololtr -i $resdir -g $masked -r $allstfile -l 80 -p -c 0.09 -s";
#say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };
#system([0..5], $find_cmd);

ok( -s $allstfile, 'Generated summary statistics for all solo-LTR matches' );
ok( -s $parsfile,  'Generated statistics for solo-LTR matches' );
ok( -s $seqfile,   'Generated sequences for all solo-LTR matches' );

my $seqct = 0;
open my $in, '<', $seqfile;
while (<$in>) { $seqct++ if /^>/; }
ok( $seqct == 14, 'Correct number of solo-LTR sequences above thresholds' );
close $in;

# clean up
remove_tree( $outdir, { safe => 1 } );
unlink $masked;

done_testing();
