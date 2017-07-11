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

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $outdir  = File::Spec->catdir($testdir, 't_family_domains');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $paired  = File::Spec->catfile($testdir, 'ref_hscan_paired.txt');
my $ext     = File::Spec->catfile($testdir, 'ref_tephra_hscan_helitrons.ext.hel.fa');
my $flank   = File::Spec->catfile($testdir, 'ref_tephra_hscan_helitrons.flanking.fa');
my $full    = File::Spec->catfile($testdir, 'ref_tephra_hscan_helitrons.hel.fa');
my $hsgff   = File::Spec->catfile($testdir, 'ref_tephra_helitrons.gff3');
my $hsfas   = File::Spec->catfile($testdir, 'ref_tephra_helitrons.fasta');
my $hslog   = File::Spec->catfile($testdir, 'ref_tephra_findhelitrons.log');

{
    my @help_args = ($cmd, 'findhelitrons', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
        #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute findhelitrons subcommand');
}

my @find_cmd = ($cmd, 'findhelitrons', '-g', $genome, '-o', $hsgff);
say STDERR join q{ }, @find_cmd;
my @ret = capture { system([0..5], @find_cmd) };
#system([0..5], @find_cmd);

my @lcvs;
find( sub { push @lcvs, $File::Find::name if -f and /.lcvs$/ }, $testdir );
#ok( @lcvs == 2, 'Correct number of termini sequences for helitron found' );
#ok( -s $paired, 'Generated paired termini sequences' );
#ok( -s $ext,    'Generated extended termini sequences' );
#ok( -s $flank,  'Generated flanking sequences' );
#ok( -s $full,   'Generated full length helitron' );

my $seqct = 0;
my ($id, $family);
open my $in, '<', $hsfas;
while (<$in>) { 
    chomp;
    if (/^>(DHH_singleton_family0)_(helitron1)_Contig57_HLAC-254L24_106214_107555/) {
	($family, $id) = ($1, $2);
	$seqct++;
    }
}
close $in;

my ($gid, $gfamily);
open my $gin, '<', $hsgff;
while (<$gin>) {
    chomp;
    next if /^#/;
    my @f = split /\t/;
    ($gid, $gfamily) = ($f[8] =~ /ID=(helitron1);family=(DHH_singleton_family0);Ontology_term=SO:0000544/);
}
close $gin;

my $tot;
open my $lin, '<', $hslog;
while (<$lin>) {
    chomp;
    if (/Results - Number of Helitron elements written.*(\d+)$/) {
	$tot = $1;
    }
}
close $lin;

ok( -e $hslog, 'findhelitrons log created' );
ok( $seqct == 1, 'Correct number of Helitrons found' );
ok( $id eq $gid, 'Correct ID added to FASTA and GFF3' );
ok( $family eq $gfamily, 'Correct family added to FASTA and GFF3' );
ok( $tot == 1, 'Correct number of elements logged' );
ok( $tot == $seqct, 'Correct number of elements logged and written');

# clean up
my @files;
find( sub { push @files, $File::Find::name if -f and /hscan/ }, $testdir );
unlink @files;
unlink $hsgff, $hsfas, $hslog;

done_testing();
