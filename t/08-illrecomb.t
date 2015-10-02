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

use Test::More tests => 8;

my $bindir = File::Spec->catdir('t', 'gt', 'bin');
local $ENV{PATH} = "$bindir:$ENV{PATH}";

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $outdir    = File::Spec->catdir($testdir, 't_family_domains');
my $resdir    = File::Spec->catdir($outdir, 'ref_ltrdigest85_combined_filtered_gypsy');
my $allstfile = File::Spec->catfile($resdir, 'gypsy_illrecomb_stats.tsv');
my $illstfile = File::Spec->catfile($resdir, 'gypsy_illrecomb_illrecstats.tsv');
my $seqfile   = File::Spec->catfile($resdir, 'gypsy_illrecomb_seqs.fasta');
my $genome    = File::Spec->catfile($testdir, 'ref.fas');
my @results   = capture { system([0..5], "$cmd illrecomb -h") };

ok(@results, 'Can execute illrecomb subcommand');

my $find_cmd = "$cmd illrecomb -i $resdir -s $allstfile -r $illstfile -o $seqfile";
say STDERR $find_cmd;

#my @ret = capture { system([0..5], $find_cmd) };
system([0..5], $find_cmd);

ok( -s $allstfile, 'Generated statistics for all gap sites' );
ok( -s $illstfile, 'Generated statistics for all putative illegetimate recombination sites' );
ok( -s $seqfile,   'Generated sequences flanking all putative illegetimate recombination sites' );

my $seqct = 0;
open my $in, '<', $seqfile;
while (<$in>) { $seqct++ if /^>/; }
ok( $seqct/2 == 16, 'Correct number of illigetimate recombination events detected' );
close $in;

my ($qmatch, $hmatch) = (0, 0);
open my $stats, '<', $illstfile;
while (<$stats>) {
    chomp;
    $qmatch++ if /^Query match string/;
    $hmatch++ if /^Hit match string/;
}

ok( $qmatch == 16, 'Correct number of illigetimate recombination events detected upstream of gap' );
ok( $hmatch == 16, 'Correct number of illigetimate recombination events detected downstream of gap' );
ok( $seqct == $qmatch+$hmatch, 'Correct number of illigetimate recombination events detected' );

unlink $allstfile, $illstfile, $seqfile;

done_testing();
