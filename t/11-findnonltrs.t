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
#use Data::Dump;

use Test::More tests => 7;

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'Ha1.fa');
my $gff     = File::Spec->catfile($testdir, 'Ha1_nonLTRs.gff3');
my $fas     = File::Spec->catfile($testdir, 'Ha1_nonLTRs.fasta');
my $outdir  = File::Spec->catdir($testdir,  'Ha1_nonLTRs');
my $log     = File::Spec->catdir($testdir,  'Ha1_tephra_findnonltrs.log');

{
    my @help_args = ($cmd, 'findnonltrs', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
        #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute findnonltrs subcommand');
}

SKIP: {
    skip 'skip lengthy tests', 6 unless $devtests;
    my @find_cmd = ($cmd, 'findnonltrs', '-g', $genome, '-o', $gff);
    say STDERR join q{ }, @find_cmd;
    my ($stdout, $stderr, @ret) = capture { system([0..5], @find_cmd) };
    #system([0..5], @find_cmd);

    ok( -e $gff, 'Can find some non-LTRs' );
    ok( -e $fas, 'Can find some non-LTRs' );

    my $seqct = 0;
    open my $in, '<', $fas;
    while (<$in>) { 
	chomp;
	if (/^>RIL/) {
	    $seqct++;
	}
    }
    close $in;

    my $gct = 0;
    open my $gin, '<', $gff;
    while (<$gin>) {
	chomp;
	next if /^#/;
	my @f = split /\t/;
	++$gct if $f[8] =~ /ID=non_LTR_retrotransposon\d+.*;family/;
    }
    close $gin;

    my $tot;
    open my $lin, '<', $log;
    while (<$lin>) {
	chomp;
	if (/Results - Number of non-LTR elements written.*\s+(\d+)$/) {
	    $tot = $1;
	}
    }
    close $lin;

    say STDERR "DEBUG: tot -> $tot";
    ok( -e $log, 'findnonltrs log created' );
    ok( $seqct == 16, 'Correct number of non-LTRs found' );
    ok( $gct == 16, 'Correct number of non-LTRs found' );
    ok( $tot == 16, 'Correct number of elements logged' );
    #ok( $tot == $seqct, 'Correct number of elements logged and written');

    ## clean up
    unlink $gff, $fas, $log;
    remove_tree( $outdir, { safe => 1 } );
}
    
done_testing();
