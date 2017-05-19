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

use Test::More tests => 5;

$| = 1;

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $gff     = File::Spec->catfile($testdir, 'ref_tirs.gff3');
my $outgff  = File::Spec->catfile($testdir, 'ref_tirs_classified.gff3');
my $outfas  = File::Spec->catfile($testdir, 'ref_tirs_classified.fasta');
my $log     = File::Spec->catfile($testdir, 'ref_tephra_classifytirs.log');

my @results = capture { system([0..5], "$cmd classifyltrs -h") };

ok(@results, 'Can execute classifytirs subcommand');

my $find_cmd = "$cmd classifytirs -g $genome -i $gff -o $outgff";
#say STDERR $find_cmd;

my @ret = capture { system([0..5], $find_cmd) };

ok( -e $outgff, 'Correctly classified TIRs' );
ok( -e $outfas, 'Correctly classified TIRs' );

my $seqct = 0;
open my $in, '<', $outfas;
while (<$in>) { $seqct++ if /^>/; }
close $in;

my $gffct = 0;
open my $gin, '<', $outgff;
while (<$gin>) { 
    chomp;
    next if /^#/;
    my @f = split /\t/;
    $gffct++ if $f[2] eq 'terminal_inverted_repeat_element'
}
close $gin;

ok( $seqct == 1, 'Correct number of TIRs classified' );
ok( $gffct == 1, 'Correct number of TIRs classified' );


## clean up
unlink $gff, $outfas, $log;
    
done_testing();
