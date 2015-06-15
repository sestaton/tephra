#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use File::Spec;

use Test::More tests => 2;

BEGIN {
    use_ok( 'Chloro' ) || print "Bail out!\n";
}

diag( "Testing Chloro $Chloro::VERSION, Perl $], $^X" );

my $cmd = File::Spec->catfile('bin', 'chloro');
ok( -x $cmd, 'Can execute chloro' );

done_testing();
