#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use HTTP::Tiny;
use HTML::TreeBuilder;
use File::Spec;
use File::Copy qw(move);
use File::Path qw(make_path);
use Tephra::Config;
use Cwd;

use Test::More tests => 10;

BEGIN {
    use_ok( 'Tephra' ) || print "Bail out!\n";
    use_ok( 'Tephra::Config' ) || print "Bail out!\n"
}

diag( "Testing Tephra $Tephra::VERSION, Perl $], $^X" );

my $cmd = File::Spec->catfile('blib', 'bin', 'tephra');
ok( -x $cmd, 'Can execute tephra' );

#chdir 't' or die $!;
#my $gt    = get_gt_exes();
#my $hscan = get_hscan();

my $cwd = getcwd();
#say STDERR $cwd and exit;
my $basedir = File::Spec->catdir($ENV{HOME}, '.tephra');
my $confobj = Tephra::Config->new( basedir => $basedir, workingdir => $cwd );
my ($gt, $hscan, $hmmbin, $moddir, $chrdir, $mgescan, $trans) = $confobj->configure_root;
my $hmmsearch = File::Spec->catfile($hmmbin, 'hmmsearch');

ok( -x $gt,        'Can execute gt for testing' );
ok( -e $hscan,     'Can execute HelitronScanner for testing' );
ok( -x $hmmsearch, 'Can execute hmmsearch' );
ok( -e $moddir,    'Configured pHMM dir for non-LTR search' );
ok( -e $chrdir,    'Configured HMM dir for non-LTR search' );
ok( -e $mgescan,   'Can build custom MGEScan for non-LTR search' );
ok( -e $trans,     'Can build translate command for non-LTR search' );

done_testing();

