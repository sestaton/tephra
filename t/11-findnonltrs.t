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

use Test::More 'no_plan';

$| = 1;

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $gff     = File::Spec->catfile($testdir, 'ref_nonLTRs.gff3');
my $fas     = File::Spec->catfile($testdir, 'ref_nonLTRs.fasta');
my $outdir  = File::Spec->catdir($testdir,  'ref_nonLTRs');
my $log     = File::Spec->catdir($testdir,  'ref_tephra_findnonltrs.log');

my $devtests = 0;
my ($exp_seqct, $exp_gct, $exp_tot) = (0, 0, 0);
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests  = 1;
    $exp_seqct = 5;
    $exp_gct   = 5;
    $exp_tot   = 5;

    $genome = File::Spec->catfile($testdir, 'TAIR10_chr1.fas');
    $gff    = File::Spec->catfile($testdir, 'TAIR10_chr1_nonLTRs.gff3');
    $fas    = File::Spec->catfile($testdir, 'TAIR10_chr1_nonLTRs.fasta');
    $outdir = File::Spec->catfile($testdir, 'TAIR10_chr1_nonLTRs');                                   
    $log    = File::Spec->catfile($testdir, 'TAIR10_chr1_tephra_findnonltrs.log');
}

{
    my @help_args = ($cmd, 'findnonltrs', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute findnonltrs subcommand');
}

my @find_cmd = ($cmd, 'findnonltrs', '-g', $genome, '-o', $gff);
say STDERR join q{ }, @find_cmd if $devtests;
my ($stdout, $stderr, @ret) = capture { system([0..5], @find_cmd) };

my ($seqct, $gct, $tot) = (0, 0, 0);
if ($devtests) {
    $exp_tot = 5;
    ok( -e $gff, 'Can find some non-LTRs' );
    ok( -e $fas, 'Can find some non-LTRs' );

    open my $in, '<', $fas;
    while (<$in>) { $seqct++ if /^>/; }
    close $in;
    say STDERR "seqct: $seqct";
    
    open my $gin, '<', $gff;
    while (<$gin>) {
	chomp;
	next if /^#/;
	my @f = split /\t/;
	++$gct if $f[8] =~ /ID=non_LTR_retrotransposon\d+.*;family/;
    }
    close $gin;
    
    open my $lin, '<', $log;
    while (<$lin>) {
	chomp;
	if (/Results - Number of non-LTR elements written.*\s+(\d+)$/) {
	    $tot = $1;
	}
    }
    close $lin;
    say STDERR "tot: $tot";
    
    ok( -e $log, 'findnonltrs log created' );
    ok( $tot == $seqct, 'Correct number of elements logged and written');
}

#say STDERR "DEBUG: $tot -> $exp_tot";
#say STDERR "DEBUG: $gct -> $exp_gct";
ok( $seqct == $exp_seqct, 'Correct number of non-LTRs found' );
ok( $gct == $exp_gct, 'Correct number of non-LTRs found' );
ok( $tot == $exp_tot, 'Correct number of elements logged' );

## clean up
unlink $gff, $fas, $log;
unlink $genome if $devtests;
unlink $genome.'.fai' if $devtests;
remove_tree( $outdir, { safe => 1 } );

#done_testing();
