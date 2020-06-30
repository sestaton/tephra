#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use Net::FTP;
use File::Spec;
use File::Copy          qw(move);
use File::Path          qw(make_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Cwd                 qw(getcwd);
use Tephra::Config::Exe;

use Test::More 'no_plan';

$| = 1;

BEGIN {
    use_ok( 'Tephra' ) || print "Bail out!\n";
    use_ok( 'Tephra::Config::Exe' ) || print "Bail out!\n"
}

diag( "Testing Tephra $Tephra::VERSION, Perl $], $^X" );

my $cmd = File::Spec->catfile('blib', 'bin', 'tephra');
ok( -x $cmd, 'Can execute tephra' );

{
    my ($stdout, $stderr, $exit) = capture { system($cmd) };
    ok($stderr, 'Can execute commands subcommand');
}

my $config = Tephra::Config::Exe->new->get_config_paths;
my ($tephrabin, $gt, $vmatch, $mkvtree, $hscan, $hmm2bin, $hmm3bin, $moddir, $chrhmm, $mgescan, $trans, $baseml, $transeq, $htslibdir, $blastn, $makeblastdb, $muscle)
    = @{$config}{qw(tephrabin gt vmatch mkvtree hscanjar hmmer2bin hmmer3bin modeldir chrhmm mgescan transcmd baseml transeq htslibdir blastn makeblastdb muscle)};

my $hmm2search = File::Spec->catfile($hmm2bin, 'hmmsearch');
my $hmm3search = File::Spec->catfile($hmm3bin, 'hmmsearch');
#my $blastn     = File::Spec->catfile($blastpath,   'blastn');
#my $vmatch     = File::Spec->catfile($tephrabin,   'vmatch');
#my $mkvtree    = File::Spec->catfile($tephrabin,   'mkvtree');

ok( -x $gt,         'Can execute gt for testing' );
ok( -x $vmatch,     'Can execute vmatch for testing' );
ok( -x $mkvtree,    'Can execute mkvtree for testing' );
ok( -x $muscle,     'Can build muscle for multi-sequence alignments' );
ok( -x $mgescan,    'Can build custom MGEScan for non-LTR search' );
ok( -x $trans,      'Can build translate command for non-LTR search' );
ok( -x $baseml,    'Can build PAML for analyzing LTR demography' );
ok( -x $transeq,    'Can build transeq for identifying coding domains' );
ok( -x $blastn,     'Can build blastn for sequence searches' );
ok( -x $makeblastdb,     'Can build blastn for sequence searches' );
ok( -x $hmm2search, 'Can execute HMMERv2 hmmsearch' );
ok( -x $hmm3search, 'Can execute HMMERv3 hmmsearch' );

ok( -e $hscan,      'Can execute HelitronScanner for testing' );
ok( -e $moddir,     'Configured pHMM dir for non-LTR search' );
ok( -e $chrhmm,     'Configured HMM dir for non-LTR search' );
ok( -e $htslibdir,  'Can build HTSlib for indexing and parsing sequence files' );

if (defined $ENV{TEPHRA_ENV} && $ENV{TEPHRA_ENV} eq 'development') {
    my $wd = getcwd();
    my $dev_file = fetch_dev_tair_data($wd);
    ok( -e $dev_file, 'Can download genome data from TAIR for testing ');
}

#done_testing();

sub fetch_dev_tair_data {
    my ($cdir) = @_;
    my $wd = File::Spec->catdir($cdir, 't', 'test_data');
    chdir $wd or die $!;

    my $host = 'ftp.arabidopsis.org';
    my $ftp = Net::FTP->new($host, Passive => 1, Debug => 0)
	or die "Cannot connect to $host: $@";

    $ftp->login or die "Cannot login ", $ftp->message;

    my $dir = '/home/tair/Sequences/whole_chromosomes/';
    $ftp->cwd($dir)
	or die "Cannot change working directory ", $ftp->message;

    my $file = 'TAIR10_chr1.fas';

    $ftp->binary();
    my $rsize = $ftp->size($file) or die "Could not get size ", $ftp->message;
    $ftp->get($file) or die "get failed ", $ftp->message;
    my $lsize = -s $file;
    die "Failed to fetch complete file: $file (local size: $lsize, remote size: $rsize)"
	unless $rsize == $lsize;

    return File::Spec->catfile($wd, $file);
}
