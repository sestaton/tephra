#!/usr/bin/env perl

use 5.010;
use strict;
use warnings FATAL => 'all';
use HTTP::Tiny;
use HTML::TreeBuilder;
use File::Spec;
use File::Copy qw(move);

use Test::More tests => 3;

BEGIN {
    use_ok( 'Tephra' ) || print "Bail out!\n";
}

diag( "Testing Tephra $Tephra::VERSION, Perl $], $^X" );

my $cmd = File::Spec->catfile('blib', 'bin', 'tephra');
ok( -x $cmd, 'Can execute tephra' );

chdir 't' or die $!;

my $host = 'http://genometools.org';
my $dir  = 'pub/binary_distributions';
my $file = 'gt_distlisting.html';
fetch_file($file, $host."/".$dir);

my $tree = HTML::TreeBuilder->new;
$tree->parse_file($file);

my ($dist, $ldist, $ldir);
for my $tag ($tree->look_down(_tag => 'a')) {
    if ($tag->attr('href')) {
	if ($tag->as_text =~ /Linux_x86_64-64bit-barebone.tar.gz\z/) {
	    $dist = $tag->as_text;
	    my $archive = join "/", $host, $dir, $dist;
	    fetch_file($dist, $archive);

	    $ldist = $dist;
	    $ldist =~ s/\.tar.gz\z//;
	    $ldir = 'gt';

	    system("tar xzf $dist") == 0 or die $!;

	    move $ldist, $ldir or die "Move failed: $!";
	    unlink $dist;
	}
    }
}

my $gt = File::Spec->catfile($ldir, 'bin', 'gt');
ok( -x $gt, 'Can execute gt for testing' );

unlink $file;

done_testing();
## methods
sub fetch_file {
    my ($file, $endpoint) = @_;
    unless (-e $file) {
	my $response = HTTP::Tiny->new->get($endpoint);
	unless ($response->{success}) {
	    die "Can't get url $endpoint -- Status: ", $response->{status}, " -- Reason: ", $response->{reason};
	}
	open my $out, '>', $file;
	print $out $response->{content};
	#sleep 1;
	close $out;
    }
}
