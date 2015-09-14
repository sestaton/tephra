package Tephra::Classify::LTRFams;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use IPC::System::Simple qw(capture);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use List::MoreUtils     qw(any none);
use Path::Class::File;
use Try::Tiny;
use Cwd;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has outdir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 1,
    coerce   => 1,
);

#has gff => (
#      is       => 'ro',
#      isa      => 'Path::Class::File',
#      required => 1,
#      coerce   => 1,
#);

#
# methods
#
sub extract_features {
    my $self = shift;
    my $fasta  = $self->genome;
    my $dir    = $self->outdir;
    my ($infile) = @_;
    
    #my ($fasta, $dir, $infile) = @_;
    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
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

    my ($start, $end, $region, $key, %feature, %ltrs);
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
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
	# full element
	my ($source, $element, $start, $end) = split /\-/, $ltrs{$ltr}{'full'};
	my $outfile = File::Spec->catfile($dir, $ltr.".fasta");
	$self->subseq($fasta, $source, $element, $start, $end, $outfile, $allfh);

	# pbs
	if ($ltrs{$ltr}{'pbs'}) {
	    my ($pbssource, $pbselement, $trna, $pbsstart, $pbsend) = split /\|\|/, $ltrs{$ltr}{'pbs'};
	    my $pbs_tmp = File::Spec->catfile($dir, $ltr."_pbs.fasta");
	    $self->subseq($fasta, $pbssource, $pbselement, $pbsstart, $pbsend, $pbs_tmp, $pbsfh);
	}

	# ppt
	if ($ltrs{$ltr}{'ppt'}) {
	    my ($pptsource, $pptelement, $pptstart, $pptend) = split /\|\|/, $ltrs{$ltr}{'ppt'};
	    my $ppt_tmp = File::Spec->catfile($dir, $ltr."_ppt.fasta");
	    $self->subseq($fasta, $source, $pptelement, $pptstart, $pptend, $ppt_tmp, $pptfh);
	}

	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltre, $s, $e) = split /\|\|/, $ltr_repeat;
	    if ($ltrct) {
		my $fiveprime_tmp = File::Spec->catfile($dir, $ltr."_5prime-ltr.fasta");
		$self->subseq($fasta, $src, $ltre, $s, $e, $fiveprime_tmp, $fivefh);
	    }
	    else {
		my $threeprime_tmp = File::Spec->catfile($dir, $ltr."_3prime-ltr.fasta");
		$self->subseq($fasta, $src, $ltre, $s, $e, $threeprime_tmp, $threfh);
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
	    $self->subseq($fasta, $src, $elem, $s, $e, $tmp, $fh);
	}
	close $fh;
	unlink $pdom_file if ! -s $pdom_file;
    }

    for my $file ($comp, $ppts, $pbs, $five_pr_ltrs, $three_pr_ltrs) {
	unlink $file if ! -s $file;
    }
}

sub subseq {
    my $self = shift;
    my ($fasta, $loc, $elem, $start, $end, $tmp, $out) = @_;
    my $cmd = "samtools faidx $fasta $loc:$start-$end > $tmp";
    $self->run_cmd($cmd);

    my $id = join "_", $loc, $elem, "$start-$end";
    if (-s $tmp) {
	my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
	while (my $seqobj = $seqio->next_seq) {
	    my $seq = $seqobj->seq;
	    if ($seq) {
		$seq =~ s/.{60}\K/\n/g;
		say $out join "\n", ">".$id, $seq;
	    }
	}
    }
    unlink $tmp;
}

__PACKAGE__->meta->make_immutable;

1;
