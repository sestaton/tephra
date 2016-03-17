package Tephra::LTR::MakeExemplars;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

=head1 NAME

Tephra::LTR::MakeExemplars - Make exemplars from a LTR retrotransposon family

=head1 VERSION

Version 0.02.1

=cut

our $VERSION = '0.02.1';
$VERSION = eval $VERSION;

has dir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 1,
    coerce   => 1,
);

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has gff => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

#
# methods
#
sub make_exemplars {
    my $self = shift;
    my $dir   = $self->dir;
    my $gff   = $self->gff;
    my $fasta = $self->genome;
    
    my $exemplars = $self->process_vmatch_args($dir);
    #dd $exemplars;

    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

    my ($start, $end, $source, $elem_id, $key, $full_feats, $exempid, 
	%feature, %ltrs, %seen, @fullseqs, @ltrseqs);
    
    for my $exemplar (keys %$exemplars) {
	my ($family) = ($exemplar =~ /(^RL[CGX]_family\d+)/);
	$exemplar =~ s/${family}_//;

	my $exemcomp = File::Spec->catfile($dir, $family."_exemplar_complete.fasta");
	my $exemltrs = File::Spec->catfile($dir, $family."_exemplar_ltrs.fasta");
	push @fullseqs, $exemcomp;
	push @ltrseqs, $exemltrs;
	
	open my $allfh, '>>', $exemcomp or die "\nERROR: Could not open file: $exemcomp\n";
	open my $ltrs_outfh, '>>', $exemltrs or die "\nERROR: Could not open file: $exemltrs\n";
	
	while (my $feature = $gffio->next_feature()) {
	    if ($feature->primary_tag eq 'LTR_retrotransposon') {
		my @string = split /\t/, $feature->gff_string;
		($elem_id) = ($string[8] =~ /ID=?\s+?(LTR_retrotransposon\d+)/);
		($source, $start, $end) = ($string[0], $feature->start, $feature->end);
		$key = join ".", $elem_id, $start, $end;
		$full_feats = join "||", $source, $feature->primary_tag, @string[3,4,6];
	    }
	    next unless defined $start && defined $end;
	    my $exemplar_id_form = join "_", $elem_id, $source;

	    if ($exemplar eq $exemplar_id_form) {
		$ltrs{$key}{'full'} = $full_feats;
		if ($feature->primary_tag eq 'long_terminal_repeat') {
		    my @string = split /\t/, $feature->gff_string;
		    if ($feature->start >= $start && $feature->end <= $end) {
			my $ltrkey = join "||", $string[0], $feature->primary_tag, @string[3,4,6];
			push @{$ltrs{$key}{'ltrs'}}, $ltrkey unless exists $seen{$ltrkey};
			$seen{$ltrkey} = 1;
		    }
		}
	    }
	}

	my %pdoms;
	my $ltrct = 0;
	for my $ltr (sort keys %ltrs) {
	    my ($element, $rstart, $rend) = split /\./, $ltr;
	    my $orient;
	    # full element
	    my ($source, $prim_tag, $start, $end, $strand) = split /\|\|/, $ltrs{$ltr}{'full'};
	    $orient = "plus" if $strand eq '+';
	    $orient = "minus" if $strand eq '-';
	    
	    my $full_tmp = File::Spec->catfile($dir, $family."_exemplar.fasta");
	    $self->subseq($fasta, $source, $element, $start, $end, $full_tmp, $allfh, $orient);

	    # ltrs
	    for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
		my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
		my $ltrs_out = File::Spec->catfile($dir, $family."_exemplar_ltrs.fasta");
		open my $ltrs_outfh, '>>', $ltrs_out;
		my $lfname = $ltr;
		my $orientation;
		if ($ltrct) {
		    $orientation = "5prime" if $strand eq '+';
		    $orientation = "3prime" if $strand eq '-';
		    $lfname .= "_$orientation-ltr.fasta" if $strand eq '+';
		    $lfname .= "_$orientation-ltr.fasta" if $strand eq '-';
		    my $fiveprime_tmp = File::Spec->catfile($dir, $lfname);
		    $self->subseq($fasta, $src, $element, $s, $e, $fiveprime_tmp, $ltrs_outfh, $orientation);
		    $ltrct = 0;
		}
		else {
		    $orientation = "3prime" if $strand eq '+';
		    $orientation = "5prime" if $strand eq '-';
		    $lfname .= "_$orientation-ltr.fasta" if $strand eq '+';
		    $lfname .= "_$orientation-ltr.fasta" if $strand eq '-';
		    my $threeprime_tmp = File::Spec->catfile($dir, $lfname);
		    $self->subseq($fasta, $src, $element, $s, $e, $threeprime_tmp, $ltrs_outfh, $orientation);
		    $ltrct++;
		}
	    }
	}
	close $allfh;
	close $ltrs_outfh;
    }

    return (\@fullseqs, @ltrseqs);
}

sub process_vmatch_args {
    my $self = shift;
    my ($dir) = @_;

    my (@fams, %exemplars);
    my $wanted  = sub { push @fams, $File::Find::name if -f and /(?:family\d+).fasta$/ };
    my $process = sub { grep ! -d, @_ };
    find({ wanted => $wanted, preprocess => $process }, $dir);

    for my $db (@fams) {
	my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
	my $index = File::Spec->catfile($path, $name."_mkvtree.index");
	my $vmerSearchOut = File::Spec->catfile($path, $name.".vmersearch");
	my $mk_args = "-db $db -dna -indexname $index -allout -pl";
	my $vm_args = "-showdesc 0 -l 20 -q $db -identity 80 $index > $vmerSearchOut";

	my $mkcmd = "mkvtree $mk_args";
	my $vmcmd = "vmatch $vm_args";
	$self->run_cmd($mkcmd);
	$self->run_cmd($vmcmd);
	my $exemplar = $self->parse_vmatch($vmerSearchOut);
	$exemplars{$exemplar} = 1;
	$self->clean_index($path);
	unlink $vmerSearchOut;
    }
    
    return \%exemplars;
}

sub parse_vmatch {
    my $self = shift;
    my ($vmerSearchOut) = @_;

    my %matches;
    open my $in, '<', $vmerSearchOut;
    while (<$in>) {
	chomp;
	next if /^#/;
	my @f = split;
	$matches{$f[1]}++;
    }
    close $in;
    
    my $tophit = (reverse sort { $matches{$a} <=> $matches{$b} } keys %matches)[0];
    return $tophit;
}

sub subseq {
    my $self = shift;
    my ($fasta, $loc, $elem, $start, $end, $tmp, $out, $orient) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($samtools) = @{$config}{qw(samtools)};
    my $cmd = "$samtools faidx $fasta $loc:$start-$end > $tmp";
    $self->run_cmd($cmd);

    my $id = join "_", $orient, $loc, $elem, $start, $end;
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

sub clean_index {
    my $self = shift;
    my ($dir) = @_;
    
    my @files;
    find( sub { push @files, $File::Find::name
		    if /\.al1|\.bck|\.bwt|\.des|\.lcp|\.llv|\.ois|\.prj|\.sds|\.skp|\.ssp|\.sti1|\.suf|\.tis/
	  }, $dir);
    unlink @files;
}

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::LTR::MakeExemplars


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
