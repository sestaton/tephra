#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Spec;
use File::Copy qw(move);

use Test::More tests => 3;

BEGIN {
    use_ok( 'Tephra' ) || print "Bail out!\n";
}

diag( "Testing Tephra $Tephra::VERSION, Perl $], $^X" );

my $cmd = File::Spec->catfile('bin', 'tephra');
ok( -x $cmd, 'Can execute tephra' );

chdir 't' or die $!;
# need to make this version agnostic
my $host = 'http://genometools.org';
my $dir  = 'pub/binary_distributions';
my $dist = 'gt-1.5.6-Linux_x86_64-64bit-barebone.tar.gz';
my $ldist = $dist;
$ldist =~ s/\.tar.gz$//;
my $ldir = 'gt';
my $archive = join "/", $host, $dir, $dist;

system("wget $archive 2>&1 > /dev/null")
    == 0 or die $!;
system("tar xzf $dist") 
    == 0 or die $!;

move $ldist, $ldir or die "Copy failed: $!";
my $gt = File::Spec->catfile($ldir, 'bin', 'gt');
ok( -x $gt, 'Can execute gt for testing' );

unlink $dist;
done_testing();
