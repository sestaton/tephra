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

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $outdir  = File::Spec->catdir($testdir, 't_family_domains');
my $resdir  = File::Spec->catdir($outdir, 'divergence_time_stats');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $gff     = File::Spec->catfile($testdir, 'ref_ltrdigest85_combined_filtered.gff3');

my @results = capture { system([0..5], "$cmd ltrage -h") };

ok(@results, 'Can execute ltrage subcommand');

my @files;
find( sub { 
    push @files, $File::Find::name if /(?:gypsy|copia|unclassified).gff3$/ }, 
      $testdir );
ok( @files == 3, 'Correctly generated subdirectories for family level classification' );

my @agefiles;
for my $file (@files) {
    (my $outfile = $file) =~ s/\.gff3/_ltrage.tsv/;
    my $find_cmd = "$cmd ltrage -g $genome -f $file -o $outfile";
    #say STDERR $find_cmd;
    my @ret = capture { system([0..5], $find_cmd) };
    push @agefiles, $outfile;
    #exit;
}

my @divfiles;
find( sub { push @divfiles, $File::Find::name if -f and /ltrage.tsv$/ }, $testdir);
ok( @divfiles == 3, 'Generated the expected number of LTR-RT age calculations');

## clean up
my @outfiles;
find( sub { push @outfiles, $File::Find::name if /^ref_ltr/ && ! /$gff/ }, $testdir);
unlink @outfiles;
unlink @agefiles;
unlink "^LTR_"; 

my @agedirs;
find( sub { push @agedirs, $File::Find::name if -d and /^ref_ltr.*ltrages$/ }, '.');
for my $dir (@agedirs) {
    remove_tree( $dir, { safe => 1 } );
}

done_testing();
