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

use Test::More tests => 3;

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $outdir  = File::Spec->catdir($testdir,  't_family_domains');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $gff     = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered_classified.gff3');
my $gindir  = File::Spec->catfile($outdir,  'ref_tephra_ltrs_combined_filtered_gypsy');
my $cindir  = File::Spec->catfile($outdir,  'ref_tephra_ltrs_combined_filtered_copia');
my $iindir  = File::Spec->catfile($outdir,  'ref_tephra_ltrs_combined_filtered_unclassified');
my @dirs = ($gindir, $cindir, $iindir);

SKIP: {
    skip 'skip development tests', 3 unless $devtests;
    {
        my @help_args = ($cmd, 'ltrage', '-h');
        my ($stdout, $stderr, $exit) = capture { system(@help_args) };
        #say STDERR "stderr: $stderr";
        ok($stderr, 'Can execute ltrage subcommand');
    }

    my $outfile = $gff;
    $outfile =~ s/\.gff3/_ltrages.tsv/;
    my @age_cmd = ($cmd, 'ltrage', '-g', $genome, '-f', $gff, '-i', $gindir, '-o', $outfile, '--all');
    my @ret = capture { system([0..5], @age_cmd) };

    ok( -s $outfile, 'Generated LTR age report for input GFF3 file');

    open my $in, '<', $outfile;
    my $h = <$in>;
    my $ct = 0;
    ++$ct while <$in>;
    close $in;
    ok( $ct == 6, 'Expected number of entries in LTR age report');

    ## clean up
    my @outfiles;
    find( sub { push @outfiles, $File::Find::name if /^ref_ltr/ && ! /$gff/ }, $testdir);
    unlink @outfiles;
};

unlink $gff;

done_testing();
