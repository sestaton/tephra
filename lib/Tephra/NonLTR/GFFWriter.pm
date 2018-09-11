package Tephra::NonLTR::GFFWriter;

use 5.014;
use Moose;
use Bio::DB::HTS::Kseq;
use File::Find;
use File::Spec;
use File::Temp qw(tempfile);
use File::Path qw(make_path remove_tree);
use Cwd        qw(abs_path);
use File::Basename;
use Sort::Naturally;
use Tephra::Annotation::Util;
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::NonLTR::GFFWriter - Take results from non-LTR search and make an annotated GFF

=head1 VERSION

Version 0.12.1

=cut

our $VERSION = '0.12.1';
$VERSION = eval $VERSION;

has genome      => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has fastadir    => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has outdir      => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has n_threshold => ( is => 'ro', isa => 'Num', required => 0, default => 0.30 );
has gff         => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );

sub write_gff {
    my $self = shift;
    my $outdir = $self->outdir;
    my $fastadir = $self->fastadir;

    my @all_clade = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1', 'R2', 'CRE');
    my $seq_dir   = File::Spec->catdir( abs_path($outdir), 'info', 'full' );
    my @cladedirs = map { File::Spec->catdir( abs_path($seq_dir), $_) } @all_clade;
    my @nonltrs;

    for my $clade (@cladedirs) {
        my $name = basename($clade);
	my $filename = $name.'.dna';
        if (-d $clade) {
	    find( sub { push @nonltrs, $File::Find::name if -f and /$filename/ }, $clade ); 	    
        }
    }

    my ($fas, $gff, $sf_elem_map) = $self->_fasta_to_gff(\@nonltrs);
    #say STDERR "Done with non-LTR search.";

    ## clean up
    my $fdir  = File::Spec->catdir( abs_path($outdir), 'f' );
    my $rdir  = File::Spec->catdir( abs_path($outdir), 'b' );
    my $rgdir = $fastadir.'_b';

    for my $dir ($fdir, $rdir, $rgdir, $fastadir) {
	remove_tree( $dir, { safe => 1 } );
    }

    return ({ fasta => $fas, gff => $gff }, $sf_elem_map);
}

sub _fasta_to_gff {
    my $self = shift;
    my ($seqs) = @_;
    my $outdir = $self->outdir;
    my $genome = $self->genome;
    my $outgff = $self->gff;

    my $util = Tephra::Annotation::Util->new;
    my $index = $self->index_ref($genome);
    my $clade_map = $self->_build_clade_map;

    my ($gname, $gpath, $gsuffix) = fileparse($outgff, qr/\.[^.]*/);
    my $tmpfname = $gname.'_tephra_nonltr_fas_XXXX';
    my $tmpgname = $gname.'_tephra_nonltr_gff_XXXX';

    my ($outf, $ffilename) = tempfile( TEMPLATE => $tmpfname, DIR => $gpath, UNLINK => 0, SUFFIX => '.fasta' );
    my ($outg, $gfilename) = tempfile( TEMPLATE => $tmpgname, DIR => $gpath, UNLINK => 0, SUFFIX => '.gff3' );

    my ($lens, $combined) = $self->_get_seq_region;

    my %regions;
    for my $file (@$seqs) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	open my $in, '<', $file or die "\n[ERROR]: Could not open file: $file\n";
	while (my $line = <$in>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		my ($id, $start, $end, $strand) = ($line =~ /^\>(.*\.?fa(?:s?t?a?))_(\d+)(?:-|_)(\d+)_(\+|\-)/);
		if (defined $id) {		    
		    $regions{$name}{$id}{$start} = join "||", $start, $end, $strand;
		}
		else {
		    say STDERR "\n[ERROR]: Could not parse sequence ID for header: '$line'. ".
			"This is a bug, please report it.\n";
		}
	    }
	}
	close $in;
    }

    my @ids = map { nsort keys %{$regions{$_}} } keys %regions;
    say $outg '##gff-version 3';

    my %seen;
    for my $seqid (@ids) {
	my ($name, $path, $suffix) = fileparse($seqid, qr/\.[^.]*/);
	if (exists $lens->{$name}) {
	    unless (exists $seen{$name}) {
		say $outg join q{ }, '##sequence-region', $name, '1', $lens->{$name};
		$seen{$name} = 1;
	    }
	}
	else {
	    say STDERR "\n[ERROR]: Could not find $name in map.\n";
	}
    }

    ##TODO: How to get the strand correct? Added in v0.03.0.
    my $ct = 0;
    my %sf_elem_map;
    for my $clade (keys %regions) {
	my $code = $clade_map->{$clade} // $clade;
	for my $seqid (nsort keys %{$regions{$clade}}) {
	    my ($name, $path, $suffix) = fileparse($seqid, qr/\.[^.]*/);
	    for my $start (sort { $a <=> $b } keys %{$regions{$clade}{$seqid}}) {
	        my ($start, $end, $strand) = split /\|\|/, $regions{$clade}{$seqid}{$start};
		my $seqname = $seqid;
		$seqname =~ s/\.fa.*//;
		#my $elem = $code."_non_LTR_retrotransposon$ct";
		my $elem = "non_LTR_retrotransposon$ct";
                my $tmp = $elem.'.fasta';
		
		my ($seq, $length) = $self->get_full_seq($index, $seqname, $start, $end);
		
		my ($filtered_seq, $adj_start, $adj_end) = $self->_filterNpercent($seq, $start, $end);
		if (defined $filtered_seq) {
		    my $id = join "_", $elem, $seqname, $adj_start, $adj_end;
		    say $outf join "\n", ">".$id, $filtered_seq;
		    say $outg join "\t", $name, 'Tephra', 'non_LTR_retrotransposon', $adj_start, $adj_end, '.', $strand, '.',
		        "ID=$elem;Ontology_term=SO:0000189";
		    my $sfcode = $util->map_superfamily_name_to_code($clade);
		    $sf_elem_map{$elem} = $sfcode;
		    $ct++;
		}
	    }
	} 
    }
    close $outf;
    close $outg;
    unlink $combined;

    return ($ffilename, $gfilename,\%sf_elem_map);
}

sub _get_seq_region {
    my $self = shift;
    my $fasdir = $self->fastadir;

    my (@seqs, %lens);
    find( sub { push @seqs, $File::Find::name if -f and /\.fa.*$/ }, $fasdir );
    my $combined = File::Spec->catfile( abs_path($fasdir), 'tephra_all_genome_seqs.fas' );
    open my $out, '>', $combined or die "\n[ERROR]: Could not open file: $combined\n";
    
    for my $seq (nsort @seqs) {
	$self->collate($seq, $out);
	my $kseq = Bio::DB::HTS::Kseq->new($seq);
	my $iter = $kseq->iterator();

	while (my $seqobj = $iter->next_seq) {
	    my $id  = $seqobj->name;
	    my $seq = $seqobj->seq;
	    my $len = length($seq);
	    $lens{$id} = $len;
	}
    }

    return (\%lens, $combined);
}

sub _filterNpercent {
    my $self = shift;
    my ($seq, $start, $end) = @_;
    my $n_thresh = $self->n_threshold;

    my ($adj_start, $adj_end);
    ## This method is for removing gap ends, which arise when going back to DNA
    ## coordinates with gappy draft genomes. Added in v0.02.7.
    my ($s) = ($seq =~ /(^N+[ATCG]{0,10}?N+?)/ig);
    my ($e) = ($seq =~ /(N+?[ATCG]{0,10}?N+$)/ig);
    if ($s) {
	my $sl = length($s);
	$adj_start = $start + $sl;
	$seq =~ s/^$s//g;
    }
    else {
	$adj_start = $start;
    }

    if ($e) {
	my $el = length($e);
	$adj_end = $end - $el;
	$seq =~ s/$e$//g;
    }
    else {
	$adj_end = $end;
    }

    my $length = length($seq); # need to get length of seq after removing gaps
    return (undef, undef, undef) unless $length;

    my $n_count = ($seq =~ tr/Nn//);
    my $n_perc  = sprintf("%.2f",$n_count/$length);

    if ($n_perc <= $n_thresh && $length >= 300) {
	$seq =~ s/.{60}\K/\n/g;
	return ($seq, $adj_start, $adj_end);
    }
    else {
	return (undef, undef, undef);
    }
}

sub _build_clade_map {
    my $self = shift;
    
    my %clade_map = (
	'CR1'    => 'RIC', # CR1 clade
	'I'      => 'RII', 
	'Jockey' => 'RIJ', 
	'L1'     => 'RIL', 
	'L2'     => 'RIL', 
	'R1'     => 'RIR', 
	'RandI'  => 'RIX', 
	'Rex'    => 'RIC', # CR1 clade, http://link.springer.com/article/10.1186/1471-2148-13-152/fulltext.html?view=classic 
	'RTE'    => 'RIT',
	'Tad1'   => 'RIX', 
	'R2'     => 'RIR', 
	'CRE'    => 'RIR', # http://www.ncbi.nlm.nih.gov/pubmed/15939396
	);

    return \%clade_map;
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::NonLTR::GFFWriter


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
