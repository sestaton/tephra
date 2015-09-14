#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Data::Printer;

my $usage  = "$0 clusterfile\n";
my $file   = shift or die $usage;
my $genome = shift or die $usage;

my $dom;
my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
if ($name =~ /(\.fa.*)/) {
    $name =~ s/$1//;
}
my (%cls, $clusnum);
open my $in, '<', $file;
while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^# args=/) {
	my ($type) = ($line =~ /\/(\S+).index\z/);
	$dom = basename($type);
	$dom =~ s/${name}_//;
	$dom =~ s/_pdom//;
	#say $dom;
	#}
	#elsif (
	    #say $line;
	#}
    }
    next if $line =~ /^# \d+/;
    if ($line =~ /^(\d+):/) {
	$clusnum = $1;
    }
    elsif ($line =~ /^\s+(\S+)/) {
	push @{$cls{$dom}{$clusnum}}, $1;
    }
}

p %cls;
