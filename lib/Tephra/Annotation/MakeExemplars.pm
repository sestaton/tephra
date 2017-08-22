package Tephra::Annotation::MakeExemplars;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use File::Find;
use File::Basename;
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use Cwd                 qw(abs_path);
use Carp 'croak';
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Logger',
     'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::Any',
     'Tephra::Role::Run::GT';

=head1 NAME

Tephra::Annotation::MakeExemplars - Make exemplars from family-level classifications

=head1 VERSION

Version 0.09.3

=cut

our $VERSION = '0.09.3';
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

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

has debug => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);
#
# methods
#
sub make_exemplars {
    my $self = shift;
    my $dir   = $self->dir->absolute->resolve;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;
   
    my $index = $self->index_ref($fasta);
    my $exemplars = $self->process_vmatch_args($dir);

    #my ($sf) = ($dir =~ /_(\w+)$/);
    my ($sf) = ($dir =~ /_((?:\w+\d+\-)?\w+)$/);
    unless (defined $sf) {
        say STDERR "\nERROR: Can't get sf from $dir $.";
    }

    my $exemcomp = File::Spec->catfile($dir, $sf.'_exemplar_complete.fasta');
    my $ltrs_out = File::Spec->catfile($dir, $sf.'_exemplar_repeats.fasta');

    open my $allfh, '>>', $exemcomp or die "\nERROR: Could not open file: $exemcomp\n";
    open my $ltrs_outfh, '>>', $ltrs_out or die "\nERROR: Could not open file: $ltrs_out\n";
    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";

    my ($source_id, $type, $strand, $exemplar_id_form, 
	$start, $end, $source, $elem_id, $key, $full_feats, 
	$exempid, %feature, %ltrs, %coord_map);

    while (my $line = <$gffio>) {
	chomp $line;
	next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	if ($feature->{type} =~ /(?:LTR|TRIM)_retrotransposon|terminal_inverted_repeat_element/) {
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

	    if ($feature->{type} =~ /long_terminal_repeat|^terminal_inverted_repeat$/) {
		my $parent = @{$feature->{attributes}{Parent}}[0];
		my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
		if ($seq_id eq $feature->{seq_id}) {
		    my ($seq_id, $type, $start, $end, $strand) =
			@{$feature}{qw(seq_id type start end strand)};
		    $strand //= '?';
		    my $ltrkey = join "||", $seq_id, $type, $start, $end, $strand, $family;
		    push @{$ltrs{$pkey}{'repeats'}}, $ltrkey;
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
	$self->subseq($index, $source, $element, $start, $end, $allfh, undef, $family);
	
	# ltrs/tirs
	for my $ltr_repeat (@{$ltrs{$ltr}{'repeats'}}) {
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
    }
    close $allfh;
    close $ltrs_outfh;

    unlink $exemcomp unless -s $exemcomp;
    unlink $ltrs_out unless -s $ltrs_out;
}

sub process_vmatch_args {
    my $self = shift;
    my ($dir) = @_;
    my $threads = $self->threads;

    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    my (@fams, %exemplars);
    my $wanted  = sub { push @fams, $File::Find::name if -f and /(?:family\d+).fasta$/ };
    my $process = sub { grep ! -d, @_ };
    find({ wanted => $wanted, preprocess => $process }, $dir);

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      my ($exemplar, $family) = @{$data_ref}{qw(exemplar family)};
			      $exemplars{$exemplar} = $family;
			} );

    for my $db (@fams) {
	$pm->start($db) and next;
	$SIG{INT} = sub { $pm->finish };
	my ($exemplar, $family) = $self->calculate_exemplars($db);
	$pm->finish(0, { exemplar => $exemplar, family => $family });
    }

    $pm->wait_all_children;
    
    return \%exemplars;
}

sub calculate_exemplars {
    my $self = shift;
    my ($db) = @_;

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $index = File::Spec->catfile($path, $name.'_mkvtree.index');
    my $vmerSearchOut = File::Spec->catfile($path, $name.'.vmersearch');
    my $mk_args = "-db $db -dna -indexname $index -allout -pl";
    my $vm_args = "-showdesc 0 -qspeedup 2 -l 20 -q $db -identity 80 $index > $vmerSearchOut";
    
    my $mkcmd = "mkvtree $mk_args";
    my $vmcmd = "vmatch $vm_args";
    say STDERR "DEBUG: $mkcmd" if $self->debug;
    say STDERR "DEBUG: $vmcmd" if $self->debug;
    $self->run_cmd($mkcmd);
    $self->run_cmd($vmcmd);
    my $exemplar = $self->parse_vmatch($vmerSearchOut);
    my ($family) = ($exemplar =~ /(^[A-Z]{3}_family\d+)/);
    $exemplar =~ s/${family}_//;
    $self->clean_index_files($index);
    unlink $vmerSearchOut;
    unlink $db;

    return ($exemplar, $family);
}

sub parse_vmatch {
    my $self = shift;
    my ($vmerSearchOut) = @_;

    my %matches;
    open my $in, '<', $vmerSearchOut or die "\nERROR: Could not open file: $vmerSearchOut\n";
    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^#/;
	$line =~ s/^\s+//;
	my @f = split /\s+/, $line;
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

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Annotation::MakeExemplars


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
