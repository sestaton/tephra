package Tephra::LTR::MakeExemplars;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use File::Find;
use File::Basename;
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use Carp 'croak';
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

=head1 NAME

Tephra::LTR::MakeExemplars - Make exemplars from a LTR retrotransposon family

=head1 VERSION

Version 0.03.4

=cut

our $VERSION = '0.03.4';
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
   
    my $index = $self->index_ref($fasta);
    my $exemplars = $self->process_vmatch_args($dir);
 
    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";

    my ($source_id, $type, $strand, $exemplar_id_form, 
	$start, $end, $source, $elem_id, $key, $full_feats, 
	$exempid, %feature, %ltrs, %coord_map);

    while (my $line = <$gffio>) {
	chomp $line;
	next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	if ($feature->{type} eq 'LTR_retrotransposon') {
	    $elem_id = @{$feature->{attributes}{ID}}[0];
	    ($source_id, $source, $type, $start, $end, $strand) 
		= @{$feature}{qw(seq_id source type start end strand)};
	    $strand //= '?';

	    $key = join "||", $elem_id, $start, $end;
	    $full_feats = join "||", $source_id, $type, $start, $end, $strand;
	    $coord_map{$elem_id} = join "||", $source_id, $start, $end;
	 
	    $exemplar_id_form = join "_", $elem_id, $source_id, $start, $end;
	}
	next unless defined $elem_id && defined $source_id && defined $start && defined $end;

	if (exists $exemplars->{$exemplar_id_form}) {
	    my $family = $exemplars->{$exemplar_id_form};
	    $ltrs{$key}{'full'} = join "||", $full_feats, $family;

	    if ($feature->{type} eq 'long_terminal_repeat') {
		my $parent = @{$feature->{attributes}{Parent}}[0];
		my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
		if ($seq_id eq $feature->{seq_id}) {
		    my ($seq_id, $type, $start, $end, $strand) =
			@{$feature}{qw(seq_id type start end strand)};
		    $strand //= '?';
		    my $ltrkey = join "||", $seq_id, $type, $start, $end, $strand, $family;
		    push @{$ltrs{$pkey}{'ltrs'}}, $ltrkey;
		}
	    }
	}
    }
    close $gffio;
 
    my %pdoms;
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
	my ($element, $rstart, $rend) = split /\|\|/, $ltr;
	my $orient;

	# full element
	my ($source, $prim_tag, $start, $end, $strand, $family) = split /\|\|/, $ltrs{$ltr}{'full'};
	my $exemcomp = File::Spec->catfile($dir, $family.'_exemplar_complete.fasta');
	my $ltrs_out = File::Spec->catfile($dir, $family.'_exemplar_ltrs.fasta');

	open my $allfh, '>>', $exemcomp or die "\nERROR: Could not open file: $exemcomp\n";
	open my $ltrs_outfh, '>>', $ltrs_out or die "\nERROR: Could not open file: $ltrs_out\n";;
	
	$self->subseq($index, $source, $element, $start, $end, $allfh, undef, $family);
	
	# ltrs
	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltrtag, $s, $e, $strand, $fam) = split /\|\|/, $ltr_repeat;
	    my $lfname = $element;
	    my $orientation;
	    if ($ltrct) {
		$orientation = '5prime' if $strand eq '+';
		$orientation = '3prime' if $strand eq '-';
		$orientation = 'unk-prime-r' if $strand eq '?';
		$self->subseq($index, $src, $element, $s, $e, $ltrs_outfh, $orientation, $family);
		$ltrct = 0;
	    }
	    else {
		$orientation = '3prime' if $strand eq '+';
		$orientation = '5prime' if $strand eq '-';
		$orientation = 'unk-prime-f' if $strand eq '?';
		$self->subseq($index, $src, $element, $s, $e, $ltrs_outfh, $orientation, $family);
		$ltrct++;
	    }
	}
	close $allfh;
	close $ltrs_outfh;
    }
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
	my ($family) = ($exemplar =~ /(^RL[CGX]_family\d+)/);
	$exemplar =~ s/${family}_//;
	$exemplars{$exemplar} = $family;
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
    my ($index, $loc, $elem, $start, $end, $out, $orient, $family) = @_;

    my $location = "$loc:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);

    croak "\nERROR: Something went wrong, this is a bug. Please report it.\n"
	unless $length;

    my $id;
    $id = join "_", $family, $elem, $loc, $start, $end if !$orient;
    $id = join "_", $orient, $family, $elem, $loc, $start, $end if $orient; # for unique IDs with clustalw

    $seq =~ s/.{60}\K/\n/g;
    say $out join "\n", ">$id", $seq;
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
