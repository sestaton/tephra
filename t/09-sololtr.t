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
use File::Copy;
#use Data::Dump;

use Test::More tests => 8;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $testfile  = File::Spec->catfile($testdir, 'RLG_family0_exemplar_ltrs.fasta');
my $outdir    = File::Spec->catdir($testdir, 't_family_domains');
my $resdir    = File::Spec->catdir($outdir, 'ref_ltrdigest85_combined_filtered_gypsy');
my $modeldir  = File::Spec->catdir($resdir, 'Tephra_LTR_exemplar_models');
my $allstfile = File::Spec->catfile($resdir, 'gypsy_sololtr_stats.tsv');
my $outfile   = File::Spec->catfile($resdir, 'gypsy_sololtrs.gff3');
my $seqfile   = File::Spec->catfile($modeldir, 
				    'RLG_family0_exemplar_ltrs_clustal-out_ref_masked99_hmmer_parsed_seq.fasta');
my $parsfile  = File::Spec->catfile($modeldir,
				    'RLG_family0_exemplar_ltrs_clustal-out_ref_masked99_hmmer_parsed.txt');
my $masked    = File::Spec->catfile($testdir, 'ref_masked99.fas');

SKIP: {
    skip 'skip development tests', 8 unless $devtests;
    copy $testfile, $resdir or die "\nERROR: copy failed $!";

    my @results   = capture { system([0..5], "$cmd sololtr -h") };
    
    ok(@results, 'Can execute sololtr subcommand');

    my $find_cmd = "$cmd sololtr -i $resdir -g $masked -r $allstfile -o $outfile -l 80 -c 0.09 -s";
    #say STDERR $find_cmd;

    my @ret = capture { system([0..5], $find_cmd) };
    #system([0..5], $find_cmd);

    ok( -s $allstfile, 'Generated summary statistics for all solo-LTR matches' );
    ok( -s $parsfile,  'Generated statistics for solo-LTR matches' );
    ok( -s $seqfile,   'Generated sequences for all solo-LTR matches' );
    
    my $seqct = 0;
    open my $in, '<', $seqfile;
    while (<$in>) { $seqct++ if /^>/; }
    ok( $seqct == 1, 'Correct number of solo-LTR sequences above thresholds' );
    close $in;

    my $soloct = 0;
    ok( -s $outfile, 'Can create GFF3 file of solo-LTRs' );
    open my $gff, '<', $outfile;
    while (<$gff>) {
	chomp;
	next if /^#/;
	my @f = split /\t/;
	$soloct++ if $f[2] eq 'solo_LTR';
    }
    close $gff;

    say STDERR "SOLOCT: $soloct";
    ok( $soloct == 1, 'Correct number of solo-LTRs found' );
    ok( $seqct == $soloct, 'Same number of sequences and elements written to GFF/FASTA' );

    # clean up
    unlink $allstfile, $outfile;
    remove_tree( $outdir, { safe => 1 } );
    unlink $masked;
};

done_testing();
