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

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $outdir  = File::Spec->catdir($testdir, 't_family_domains');
my $resdir  = File::Spec->catdir($outdir, 'divergence_time_stats');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $gff     = File::Spec->catfile($outdir, 'ref_ltrdigest85_combined_filtered_families.gff3');
my $gindir  = File::Spec->catfile($outdir, 'ref_ltrdigest85_combined_filtered_gypsy');
my $cindir  = File::Spec->catfile($outdir, 'ref_ltrdigest85_combined_filtered_copia');
my $iindir  = File::Spec->catfile($outdir, 'ref_ltrdigest85_combined_filtered_unclassified');
my @dirs = ($gindir, $cindir, $iindir);

SKIP: {
    skip 'skip development tests', 2 unless $devtests;
    my @results = capture { system([0..5], "$cmd ltrage -h") };

    ok(@results, 'Can execute ltrage subcommand');
    
    my $outfile = $gff;
    $outfile =~ s/\.gff3/_ltrages.tsv/;
    my $age_cmd = "$cmd ltrage -g $genome -f $gff -i $gindir -o $outfile";
    #say STDERR $age_cmd;
    my @ret = capture { system([0..5], $age_cmd) };

    ok( -s $outfile, 'Generated LTR age report for input GFF3 file');

    my @resdirs;
    find( sub { push @resdirs, $File::Find::name if -d and /ltrages/ }, $testdir);
    ok( @resdirs == 1, 'Generated the expected number of LTR-RT age calculations');
    
    ## clean up
    my @outfiles;
    find( sub { push @outfiles, $File::Find::name if /^ref_ltr/ && ! /$gff/ }, $testdir);
    unlink @outfiles;
    
    #my @agedirs;
    #find( sub { push @agedirs, $File::Find::name if -d and /^ref_ltr.*ltrages$/ }, '.');
    for my $dir (@resdirs) {
	remove_tree( $dir, { safe => 1 } );
    }
};

done_testing();
