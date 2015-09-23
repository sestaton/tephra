#!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Data::Printer;
use Sort::Naturally;
use List::MoreUtils qw(indexes any natatime);
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
	my $element = $1;
	$element =~ s/_\d+-\d+$//;
	push @{$cls{$dom}{$clusnum}}, $element;
    }
}

my (%elem_sorted, %multi_cluster_elem);
for my $pdom (keys %cls) {
    for my $clsnum (keys %{$cls{$pdom}}) {
	for my $elem (@{$cls{$pdom}{$clsnum}}) {
	    push @{$elem_sorted{$elem}}, { $pdom => $clsnum };
	}
    }
}

#print join q{ }, (nsort keys %cls);
#print "\n";

my %dom_orgs;
for my $element (keys %elem_sorted) {
    my $string;
    my %pdomh;
    @pdomh{keys %$_} = values %$_ for @{$elem_sorted{$element}};
    for my $pdom (nsort keys %cls) {
	if (exists $pdomh{$pdom}) {
	    $string .= length($string) ? "|$pdomh{$pdom}" : $pdomh{$pdom};
	}
	else {
	    $string .= length($string) ? "|N" : "N";
	}
    }
    push @{$dom_orgs{$string}}, $element;
    undef $string;
}

dd \%dom_orgs and exit;

my (%ind, %idx_id);
my %joins;

for my $org (sort keys %dom_orgs) {
    my @ar = split /\|/, $org;
    my @i  = indexes { $_ =~ /\d+/ } @ar;  # indexes with a domain cluster id
    my $k  = join "||", @i;                # as a string
    my $vc = @i;                           # the count of domains assigned to a cluster
    my $v  = join "||", @ar[@i];           # the cluster IDs as a string

    if (@i > 1) {
	push @{$ind{$k}}, { $v => $org };
	push @{$idx_id{$k}}, $org
    }
}

#my %seen;
#my $it = natatime 2, (keys %ind);
#dd \(keys %ind);
my $ct = 0;
for my $i (keys %ind) { 
    for my $j (keys %ind) {
	next if $i eq $j;
	my @i_arr = split /\|\|/, $i;
	my @j_arr = split /\|\|/, $j;

	for my $clsidi (@i_arr) {
	    for my $clsidj (@j_arr) {
		$ct++ if $clsidi == $clsidj;
	    }
	}
	
	if ($ct > 1) {
	    my @ids;
	    #my @cls_id_nums = ($i =~ /(\d+)/g);
	    if (exists $joins{$i} && exists $idx_id{$i}) {
		#for my $num (@cls_id_nums) {
		    #push @ids, $_ if any { $num == $_ } 
		    #map { /\d+/; split /\|/ } @{$idx_id{$i}};
		#}
		push @{$joins{$i}}, $j; #@{$idx_id{@ids}};
	    }
	    elsif (exists $joins{$j} && exists $idx_id{$j}) {
		#for my $num (@cls_id_nums) {
		    #push @ids, $_ if any { $num == $_ } 
		    #map { /\d+/; split /\|/ } @{$idx_id{$j}};
		#}		
		push @{$joins{$j}}, $i; #@{$idx_id{@ids}};
	    }
	    else {
		#for my $num (@cls_id_nums) {
		    #push @ids, $_ if any { $num == $_ } 
		    #map { /\d+/; split /\|/ } @{$idx_id{$i}};
		#}		
		push @{$joins{$i}}, $j; #@{$idx_id{@ids}};
	    }
	}
	#@cls_id_nums = ();
	#@ids = ();
	$ct  = 0;
    }
}

dd \%ind;
dd \%idx_id;
dd \%joins;

exit;
