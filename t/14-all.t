#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use autodie             qw(open);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use File::Path          qw(make_path remove_tree);
use File::Spec;
use File::Copy;
use File::Find;

use Test::More tests => 4;

$| = 1;

my $devtests = 0;
if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    $devtests = 1;
}

my $cmd       = File::Spec->catfile('blib', 'bin', 'tephra');
my $testdir   = File::Spec->catdir('t', 'test_data');
my $genome    = File::Spec->catfile($testdir, 'ref.fas');
my $thrlog    = File::Spec->catfile($testdir, 'tephra_fragment_searches.log');
my $log       = File::Spec->catfile($testdir, 'tephra_full.log');
my $gff       = File::Spec->catfile($testdir, 'tephra_transposons.gff3');
my $fas       = File::Spec->catfile($testdir, 'tephra_transposons.fasta');
my $finaldb   = File::Spec->catfile($testdir, 'tephra_transposons_complete.fasta');
#my $ctestfile = File::Spec->catfile($testdir, 'tephra_copia_exemplar_repeats.fasta');
#my $gtestfile = File::Spec->catfile($testdir, 'tephra_gypsy_exemplar_repeats.fasta');
#my $ltrcdir   = File::Spec->catdir($testdir,  'ref_tephra_ltrs_trims_classified_results');
#my $cresdir   = File::Spec->catdir($ltrcdir,  'ref_tephra_ltrs_copia');
#my $gresdir   = File::Spec->catdir($ltrcdir,  'ref_tephra_ltrs_gypsy');

{
    my @help_args = ($cmd, 'all', '-h');
    my ($stdout, $stderr, $exit) = capture { system(@help_args) };
    #say STDERR "stderr: $stderr";
    ok($stderr, 'Can execute all subcommand');
}

SKIP: {
    skip 'skip lengthy tests', 3 unless $devtests;

    my $config = write_config($testdir);
    #make_path( $ltrcdir, {verbose => 0, mode => 0771,} );
    #make_path( $cresdir, {verbose => 0, mode => 0771,} );
    #make_path( $gresdir, {verbose => 0, mode => 0771,} );
    #copy $ctestfile, $cresdir or die "\n[ERROR]: copy failed $!";
    #copy $gtestfile, $gresdir or die "\n[ERROR]: copy failed $!";
    
    my @all_cmd = ($cmd, 'all', '-c', $config);
    say STDERR join q{ }, @all_cmd;
    #my ($astdout, $astderr, @aret) = capture { system([0..5], $all_cmd) };
    system([0..5], @all_cmd);

    ok( -e $gff, 'Can run full tephra pipeline and generate combined GFF3' );
    ok( -e $fas, 'Can run full tephra pipeline and generate combined FASTA' );
    ok( -e $log, 'Can run full tephra pipeline and log results' );

    ## clean up
    unlink $gff, $fas, $log, $thrlog, $config, $finaldb;

    my @outfiles;
    find( sub { push @outfiles, $File::Find::name if /^ref_|\.fai$/ }, $testdir);
    for my $res (@outfiles) {
	unlink $res 
	    if -f $res;
	remove_tree( $res, { safe => 1 } ) 
	    if -d $res;
    }
}
unlink $genome.'.fai' 
    if -e $genome.'.fai';

system('gunzip', "$testdir/ref.fas.gz") == 0 or die $!;

done_testing();

#
# methods
#
sub write_config {
    my ($testdir) = @_;
    my $config = File::Spec->catfile($testdir, 'tephra_all_config.yml');

    system('gzip', "$testdir/ref.fas") == 0 or die $!;

    my $conf = "all:
  - logfile:          $testdir/tephra_full.log
  - genome:           $testdir/ref.fas.gz
  - outfile:          $testdir/tephra_transposons.gff3
  - repeatdb:         $testdir/repdb.fas 
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
  - all:              YES
maskref:
  - percentid:        50
  - hitlength:        50
  - splitsize:        5000
  - overlap:          100
sololtr:
  - percentid:        39
  - percentcov:       80
  - matchlen:         80
  - numfamilies:      20
  - allfamilies:      NO
tirage:
  - all:              YES";

    open my $out, '>', $config;
    say $out $conf;
    close $out;
    
    return $config;
}
