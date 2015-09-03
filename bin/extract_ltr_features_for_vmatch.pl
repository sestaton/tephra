#!/usr/bin/env perl

## features: LTRs, PBS/PPT, protein_domains
##TODO: need to add feature name, element source, to header

use 5.010;
use strict;
use warnings;
use autodie qw(open);
use File::Basename;
use File::Find;
use File::Path qw(make_path);
use IPC::System::Simple qw(system);
use Try::Tiny;
use Getopt::Long;
use Bio::Tools::GFF;
use Cwd;
use Data::Dump;
use experimental 'signatures';

my %opt;
my %ltrs;
my $ltrct = 0;
my $hash;

GetOptions(\%opt, 'infile|i=s', 'fasta|f=s', 'dir|d=s');

usage() and exit(0) if !$opt{infile} or !$opt{fasta} or !$opt{dir};
die "\nERROR: '$opt{dir}' already exists. Exiting." if -d $opt{dir};

unless ( -d $opt{dir} ) {
    make_path( $opt{dir}, {verbose => 0, mode => 0771,} );
}

extract_features($opt{fasta}, $opt{dir}, $opt{infile});

sub extract_features ($fasta, $dir, $infile) {
    my ($name, $path, $suffix) = fileparse($fasta, qr/\.[^.]*/);
    my $comp = File::Spec->catfile($dir, $name."_complete.fasta");
    my $ppts = File::Spec->catfile($dir, $name."_ppt.fasta");
    my $pbs  = File::Spec->catfile($dir, $name."_pbs.fasta");
    my $five_pr_ltrs  = File::Spec->catfile($dir, $name."_5prime-ltrs.fasta");
    my $three_pr_ltrs = File::Spec->catfile($dir, $name."_3prime-ltrs.fasta");

    open my $allfh, '>>', $comp;
    open my $pptfh, '>>', $ppts;
    open my $pbsfh, '>>', $pbs;
    open my $fivefh, '>>', $five_pr_ltrs;
    open my $threfh, '>>', $three_pr_ltrs;
    
    my $gffio = Bio::Tools::GFF->new( -file => $infile, -gff_version => 3 );

    my ($start, $end, $region, $key, %feature);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'LTR_retrotransposon') {
	    my @string = split /\t/, $feature->gff_string;
	    ($region) = ($string[8] =~ /ID=?\s+?(LTR_retrotransposon\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $key = join ".", $region, $start, $end;
	    $ltrs{$key}{'full'} = join "-", $string[0], $feature->primary_tag, @string[3..4];
	}
	next unless defined $start && defined $end;
	if ($feature->primary_tag eq 'long_terminal_repeat') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		push @{$ltrs{$key}{'ltrs'}}, 
	            join "||", $string[0], $feature->primary_tag, @string[3..4];
	    }
	}
	elsif ($feature->primary_tag eq 'primer_binding_site') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my ($name) = ($string[8] =~ /trna \"?(\w+.*)\"?\s+\;/);
		$ltrs{$key}{'pbs'} = 
		    join "||", $string[0], $feature->primary_tag, $name, @string[3..4];
	    }
	}
	elsif ($feature->primary_tag eq 'protein_match') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my ($name) = ($string[8] =~ /name \"?(\w+)\"?/);
		push @{$ltrs{$key}{'pdoms'}},
		    join "||", $string[0], $feature->primary_tag, $name, @string[3..4];
	    }
	}
	elsif ($feature->primary_tag eq 'RR_tract') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		$ltrs{$key}{'ppt'} =
		    join "||", $string[0], $feature->primary_tag, @string[3..4];
	    }
	}
    }

    my %pdoms;
    
    for my $ltr (sort keys %ltrs) {
	# full element
	my ($source, $element, $start, $end) = split /\-/, $ltrs{$ltr}{'full'};
	my $outfile = File::Spec->catfile($dir, $ltr.".fasta");
	subseq($fasta, $source, $element, $start, $end, $outfile, $allfh);

	# pbs
	if ($ltrs{$ltr}{'pbs'}) {
	    my ($pbssource, $pbselement, $trna, $pbsstart, $pbsend) = split /\|\|/, $ltrs{$ltr}{'pbs'};
	    my $pbs_tmp = File::Spec->catfile($dir, $ltr."_pbs.fasta");
	    subseq($fasta, $pbssource, $pbselement, $pbsstart, $pbsend, $pbs_tmp, $pbsfh);
	}
	
	# ppt
	if ($ltrs{$ltr}{'ppt'}) {
	    my ($pptsource, $pptelement, $pptstart, $pptend) = split /\|\|/, $ltrs{$ltr}{'ppt'};
	    my $ppt_tmp = File::Spec->catfile($dir, $ltr."_ppt.fasta");
	    subseq($fasta, $source, $pptelement, $pptstart, $pptend, $ppt_tmp, $pptfh);
	}
	
	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltre, $s, $e) = split /\|\|/, $ltr_repeat;
	    if ($ltrct) { 
		my $fiveprime_tmp = File::Spec->catfile($dir, $ltr."_5prime-ltr.fasta");
		subseq($fasta, $src, $ltre, $s, $e, $fiveprime_tmp, $fivefh);
	    }
	    else {
		my $threeprime_tmp = File::Spec->catfile($dir, $ltr."_3prime-ltr.fasta");
		subseq($fasta, $src, $ltre, $s, $e, $threeprime_tmp, $threfh);
		$ltrct++;
	    }
	}
	$ltrct = 0;
	
	if ($ltrs{$ltr}{'pdoms'}) {
	    for my $ltr_repeat (@{$ltrs{$ltr}{'pdoms'}}) { 
		my ($src, $what, $name, $s, $e ) = split /\|\|/, $ltr_repeat;
		#"Ha10||protein_match||UBN2||132013916||132014240",
		push @{$pdoms{$name}}, join "||", $src, $element, $s, $e;
	    }
	}
    }
    close $allfh;
    close $pptfh;
    close $pbsfh;
    close $fivefh;
    close $threfh;
    
    for my $pdom_type (keys %pdoms) {
	my $pdom_file = File::Spec->catfile($dir, $pdom_type."_pdom.fasta");
	open my $fh, '>>', $pdom_file;
	for my $ltrpdom (@{$pdoms{$pdom_type}}) {
	    my ($src, $elem, $s, $e) = split /\|\|/, $ltrpdom;
	    my $tmp = File::Spec->catfile($dir, $elem."_".$pdom_type.".fasta");
	    subseq($fasta, $src, $elem, $s, $e, $tmp, $fh);
	}
	close $fh;
    }

    
}

sub subseq ($fasta, $loc, $elem, $start, $end, $tmp, $out) {
    my $cmd = "samtools faidx $fasta $loc:$start-$end > $tmp";
    try {
	system([0..5], $cmd);
    }
    catch {
	die "\nERROR: $cmd failed on $tmp. Here is the exception: $_\n";
    };

    my $id = join "_", $loc, $elem, "$start-$end";
    my @aux = undef;
    my ($name, $comm, $seq, $qual);
    open my $in, '<', $tmp;
    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	$seq =~ s/.{60}\K/\n/g;
	say $out join "\n", ">".$id, $seq;
    }
    close $in;
    unlink $tmp;
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}

sub usage {
    my $script = basename($0);
    print STDERR <<END
	
USAGE: $script -i file.gff -f seqs.fas -d dirname

Required:
 -i|infile    :    GFF file to extract gene coordinates from
 -f|fasta     :    FASTA file to pull the gene regions from.
 -d|dir       :    A directory name to place the resulting FASTA files.
    
Options:
 -h|help      :    Print usage statement (not implemented).
 -m|man       :    Print full documentation (not implemented).

END
}
