#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Data::Printer;
use Sort::Naturally;
use List::MoreUtils qw(indexes any);
use Data::Dump;

my $usage  = "$0 clusterfile\n";
my $file   = shift or die $usage;
my $genome = shift or die $usage;

my $dom;
my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
if ($name =~ /(\.fa.*)/) {
    $name =~ s/$1//;
}
my (%cls, %all_seqs, %all_pdoms, $clusnum);
open my $in, '<', $file;
while (my $line = <$in>) {
    chomp $line;
    if ($line =~ /^# args=/) {
	my ($type) = ($line =~ /\/(\S+).index\z/);
	$dom = basename($type);
	$dom =~ s/${name}_//;
	$dom =~ s/_pdom//;
	$all_pdoms{$dom} = 1;
    }
    next if $line =~ /^# \d+/;
    if ($line =~ /^(\d+):/) {
	$clusnum = $1;
    }
    elsif ($line =~ /^\s+(\S+)/) {
	push @{$cls{$dom}{$clusnum}}, $1;
    }
}

#dd \%cls and exit;

my (%elem_sorted, %multi_cluster_elem);
for my $pdom (keys %cls) {
    for my $clsnum (keys %{$cls{$pdom}}) {
	for my $elem (@{$cls{$pdom}{$clsnum}}) {
	    push @{$elem_sorted{$elem}}, { $pdom => $clsnum };
	}
    }
}

#dd %elem_sorted and exit;

#print join q{ }, (nsort keys %cls);
#print "\n";

my %dom_orgs;
for my $element (keys %elem_sorted) {
    my $string;
    for my $pdomh (@{$elem_sorted{$element}}) {
	for my $pdom (nsort keys %cls) {
	    if (exists $pdomh->{$pdom}) {
		$string .= length($string) ? "|$pdomh->{$pdom}" : $pdomh->{$pdom};
	    }
	    else {
		$string .= length($string) ? "|N" : "N";
	    }
	}
	push @{$dom_orgs{$string}}, $element;
	undef $string;
    }
}

#dd %dom_orgs and exit;

say join "\t", "DomainOrg", "DomainCt", "IndexOffsets", "ClusterIDs(values)";
my %ind;
for my $org (sort keys %dom_orgs) {
    #say $org;
    my @ar = split /\|/, $org;
    my @i  = indexes { $_ =~ /\d+/ } @ar;  # indexes with a domain cluster id
    my $k  = join "||", @i;                # as a string
    my $vc = @i;                           # the count of domains assigned to a cluster
    my $v  = join "||", @ar[@i];           # the cluster IDs as a string

    if (exists $ind{$k}) {
	#say "same position: $k";
	#say $org;
	#say join "\n", @{$ind{$k}};
	#say "same cluster: $v";
    }
    else {
	push @{$ind{$k}}, $org;
    }
}

#dd \%ind;

exit;
