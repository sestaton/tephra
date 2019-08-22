#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open system);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use List::Util          qw(sum);
use File::Find;
use File::Spec;
#use Data::Dump::Color;

use Test::More tests => 19;

$| = 1;

my $cmd     = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir = File::Spec->catdir('t', 'test_data');
my $genome  = File::Spec->catfile($testdir, 'ref.fas');
my $outgff  = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered.gff3');
my $outfas  = File::Spec->catfile($testdir, 'ref_tephra_ltrs_combined_filtered.fasta');
my $log     = File::Spec->catfile($testdir, 'ref_tephra_findltrs.log');

## these are subsets for testing
#my $model   = File::Spec->catfile($testdir, 'te.hmm');
#my $trnas   = File::Spec->catfile($testdir, 'trnas.fas');

my $config  = write_config($testdir, $log);
ok( -e $config, 'Can create config file for testing' );

{
    my @help_args = ($cmd, 'findltrs', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute findltrs subcommand');
}

my @find_args = ($cmd, 'findltrs', '-c', $config); # == 0 or die $!;
#say STDERR join q{ }, @find_args;

my ($stdout, $stderr, $exit) = capture { system(@find_args) }; 
ok( -e $outgff, 'Can find some LTRs' );

for my $line (split /^/, $stderr) {
    if ($line =~ /combined elements discovered \(prior to filtering steps\)/) {
        my ($tot) = $line =~ /(\d+)$/;
        ok( $tot == 7, 'Correct number of elements discovered prior to filtering' );
    }
    elsif ($line =~ /removed by matches to gene set/) {
        my ($rm_filt) = $line =~ /(\d+)$/;
        ok( $rm_filt == 0, 'Correct number of elements filtered by matches to gene set' );
    }
    elsif ($line =~ /input to additional filtering steps/) {
        my ($in_filt) = $line =~ /(\d+)$/;
        ok( $in_filt == 7, 'Correct number of elements input to additional filtering steps' );
    }
    elsif ($line =~ /length_filtered/) {
	my ($l_filt) = $line =~ /(\d+)$/;
	#say STDERR "l_filt: $l_filt";
	ok( $l_filt == 0, 'Correct number of elements filtered by length' );
    }
    elsif ($line =~ /compound_gyp_cop_filtered/) {
	my ($gc_filt) = $line =~ /(\d+)$/;
	#say STDERR "gc_filt: $gc_filt";
	ok( $gc_filt == 0, 'Correct number of elements filtered compound gypsy/copia' );
    }
    elsif ($line =~ /n_perc_filtered/) {
	my ($n_filt) = $line =~ /(\d+)$/;
	#say STDERR "n_filt: $n_filt";
	ok( $n_filt == 0, 'Correct number of elements filtered by N-percentage' );
    }
    elsif ($line =~ /\'relaxed\' constraints/) {
	my ($r_ct) = $line =~ /(\d+)$/;
	#say STDERR "r_ct: $r_ct";
	ok( $r_ct == 5, 'Correct number of combined elements found by relaxed constraints' );
    }
    elsif ($line =~ /\'strict\' constraints/) {
	my ($s_ct) = $line =~ /(\d+)$/;
	#say STDERR "s_ct: $s_ct";
	ok( $s_ct == 2, 'Correct number of total elements found by strict constraints' );
    }
    elsif ($line =~ /\'best\'/) {
	my ($b_ct) = $line =~ /(\d+)$/;
	#say STDERR "b_ct: $b_ct";
	ok( $b_ct == 1, 'Correct number of best elements found' );
    }
    elsif ($line =~ /\'combined\'/) {
	my ($c_ct) = $line =~ /(\d+)$/;
	#say STDERR "c_ct: $c_ct";
	ok( $c_ct == 6, 'Correct number of combined elements found' );
    }
    elsif ($line =~ /Total elements written/) {
	my ($t_ct) = $line =~ /(\d+)$/;
	#say STDERR "t_ct: $t_ct";
	ok( $t_ct == 6, 'Correct number of total elements found' );
    }
}

ok( -e $outgff, 'Correctly classified LTRs' );
ok( -e $outfas, 'Correctly classified LTRs' );
    
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
    $gffct++ if $f[2] eq 'LTR_retrotransposon'
}
close $gin;

ok( $seqct == 6, 'Correct number of LTRs classified' );
ok( $gffct == 6, 'Correct number of LTRs classified' );
ok( -e $log,     'Logfile of results generated'      );

## clean up
my @outfiles;
find( sub { push @outfiles, $File::Find::name if /^ref_/ && ! /combined_filtered.gff3/ }, $testdir);
unlink @outfiles;
unlink $config, $log;
    
done_testing();

sub write_config {
    my ($testdir, $log) = @_;
    my $config = File::Spec->catfile($testdir, 'tephra_ltr_config.yml');

    my $conf = "all:
  - logfile:          $log
  - genome:           $testdir/ref.fas
  - outfile:          $testdir/tephra_transposons.gff3
  - repeatdb:         $testdir/repdb.fas 
  - genefile:         $testdir/devtest_gene_seqs.fas
  - trnadb:           TephraDB
  - hmmdb:            TephraDB
  - threads:          2
  - clean:            YES
  - debug:            NO
  - subs_rate:        1e-8
findltrs:
  - dedup:            NO
  - tnpfilter:        NO
  - domains_required: NO
  - ltrharvest:
    - mintsd: 4
    - maxtsd: 6
    - minlenltr: 100
    - maxlenltr: 6000
    - mindistltr: 1500
    - maxdistltr: 25000
    - seedlength: 30
    - tsdradius: 60
    - xdrop: 5
    - swmat: 2 
    - swmis: -2
    - swins: -3
    - swdel: -3
    - overlaps: best
  - ltrdigest:
    - pptradius: 30
    - pptlen: 8 30
    - pptagpr: 0.25
    - uboxlen: 3 30
    - uboxutpr: 0.91
    - pbsradius: 30
    - pbslen: 11 30
    - pbsoffset: 0 5
    - pbstrnaoffset: 0 5
    - pbsmaxeditdist: 1
    - pdomevalue: 10E-6
    - pdomcutoff: NONE
    - maxgaplen: 50
classifyltrs:
  - percentcov:       50
  - percentid:        80
  - hitlen:           80
illrecomb:
  - repeat_pid:       10
ltrage:
  - all:              NO
maskref:
  - percentid:        80
  - hitlength:        70
  - splitsize:        5000000
  - overlap:          100
sololtr:
  - percentid:        39
  - percentcov:       80
  - matchlen:         80
  - numfamilies:      20
  - allfamilies:      NO
tirage:
  - all:              NO";

    open my $out, '>', $config;
    say $out $conf;
    close $out;
    
    return $config;
}
