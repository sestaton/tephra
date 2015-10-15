#!/usr/bin/env perl

use 5.010;
use strict;
use warnings;
use Bio::SeqIO;
use Bio::SeqUtils;
use Data::Dump;

my $in = shift or die $!;
my $out = shift or die $!;

my $seqin  = Bio::SeqIO->new(-file => $in, -format => 'fasta');
my $seqout = Bio::SeqIO->new(-file => ">>$out", -format => 'fasta'); 

my $seqobj = $seqin->next_seq;
for my $frame (0..2) { ## forward 3 frames
    my $prot_obj = $seqobj->translate(-frame => $frame);
    my $id = $prot_obj->id;
    $id .= "_$frame";
    $prot_obj->id($id);
    $seqout->write_seq($prot_obj);
}

sub get_6frames {
    my ($seqin, $seqout) = @_;
    my $seqobj = $seqin->next_seq;
    my @seqs = Bio::SeqUtils->translate_6frames($seqobj); 
    for my $psobj (@seqs) {
	$seqout->write_seq($psobj);
    }
}
