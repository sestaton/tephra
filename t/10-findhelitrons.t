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

use Test::More tests => 7;

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $outdir    = File::Spec->catdir($testdir, 't_family_domains');
my $genome    = File::Spec->catfile($testdir, 'ref.fas');
my $paired    = File::Spec->catfile($testdir, 'ref_hscan_paired.txt');
my $ext       = File::Spec->catfile($testdir, 'ref_tephra_hscan_helitrons.ext.hel.fa');
my $flank     = File::Spec->catfile($testdir, 'ref_tephra_hscan_helitrons.flanking.fa');
my $full      = File::Spec->catfile($testdir, 'ref_tephra_hscan_helitrons.hel.fa');
my $hsgff     = File::Spec->catfile($testdir, 'ref_tephra_helitrons.gff3');

my @results   = capture { system([0..5], "$cmd findhelitrons -h") };
ok(@results, 'Can execute findhelitrons subcommand');

my $find_cmd = "$cmd findhelitrons -g $genome -o $hsgff";
#say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };
#system([0..5], $find_cmd);

my @lcvs;
find( sub { push @lcvs, $File::Find::name if -f and /.lcvs$/ }, $testdir );
ok( @lcvs == 2, 'Correct number of termini sequences for helitron found' );
ok( -s $paired, 'Generated paired termini sequences' );
ok( -s $ext,    'Generated extended termini sequences' );
ok( -s $flank,  'Generated flanking sequences' );
ok( -s $full,   'Generated full length helitron' );

my $seqct = 0;
open my $in, '<', $full;
while (<$in>) { $seqct++ if /^>/; }
ok( $seqct == 1, 'Correct number of Helitrons found' );
close $in;

# clean up
my @files;
find( sub { push @files, $File::Find::name if -f and /hscan/ }, $testdir );
unlink @files;
unlink $hsgff;

done_testing();
