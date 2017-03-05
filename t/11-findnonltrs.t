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

use Test::More tests => 3;

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

my @results = capture { system([0..5], "$cmd findnonltrs -h") };
ok( @results, 'Can execute findnonltrs subcommand' );

SKIP: {
    skip 'skip lengthy tests', 2 unless $devtests;
    my $find_cmd = "$cmd findnonltrs -g $genome -o $gff";
    #say STDERR $find_cmd;

    #my ($stdout, $stderr, @ret) = capture { system([0..5], $find_cmd) };
    system([0..5], $find_cmd);

    ok( -e $gff, 'Can find some non-LTRs' );
    ok( -e $fas, 'Can find some non-LTRs' );

    ## clean up
    unlink $gff, $fas;
    remove_tree( $outdir, { safe => 1 } );
}
    
done_testing();
