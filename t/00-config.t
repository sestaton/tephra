#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use HTTP::Tiny;
use HTML::TreeBuilder;
use File::Spec;
use File::Copy qw(move);
use File::Path qw(make_path);
use Tephra::Config::Exe;

use Test::More tests => 15;

BEGIN {
    use_ok( 'Tephra' ) || print "Bail out!\n";
    use_ok( 'Tephra::Config::Exe' ) || print "Bail out!\n"
}

diag( "Testing Tephra $Tephra::VERSION, Perl $], $^X" );

my $cmd = File::Spec->catfile('blib', 'bin', 'tephra');
ok( -x $cmd, 'Can execute tephra' );

my $config = Tephra::Config::Exe->new->get_config_paths;
my ($gt, $hscan, $hmm2bin, $hmm3bin, $moddir, $chrdir, $mgescan, $trans, $pamlbin, $transeq, $htslibdir, $blast)
    = @{$config}{qw(gt hscanjar hmmer2bin hmmer3bin modeldir hmmdir mgescan transcmd pamlbin transeq htslibdir blastpath)};

my $hmm2search = File::Spec->catfile($hmm2bin, 'hmmsearch');
my $hmm3search = File::Spec->catfile($hmm3bin, 'hmmsearch');
my $blastn     = File::Spec->catfile($blast,  'blastn');

ok( -x $gt,         'Can execute gt for testing' );
ok( -e $hscan,      'Can execute HelitronScanner for testing' );
ok( -x $hmm2search, 'Can execute HMMERv2 hmmsearch' );
ok( -x $hmm3search, 'Can execute HMMERv3 hmmsearch' );
ok( -e $moddir,     'Configured pHMM dir for non-LTR search' );
ok( -e $chrdir,     'Configured HMM dir for non-LTR search' );
ok( -e $mgescan,    'Can build custom MGEScan for non-LTR search' );
ok( -e $trans,      'Can build translate command for non-LTR search' );
ok( -e $pamlbin,    'Can build paml for analyzing LTR demography' );
ok( -e $transeq,    'Can build transeq for identify coding domains' );
ok( -e $blastn,     'Can build blastn for sequence searches' );
ok( -e $htslibdir,  'Can build HTSlib for indexing and parsing sequence files' );

done_testing();

