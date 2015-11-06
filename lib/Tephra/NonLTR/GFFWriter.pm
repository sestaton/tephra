package Tephra::NonLTR::GFFWriter;

use 5.010;
use Moose;
use Bio::SeqIO;
use File::Find;
use File::Spec;
use File::Path qw(make_path);
use File::Basename;
use Sort::Naturally;
use Data::Printer;
use namespace::autoclean;

=head1 NAME

Tephra::NonLTR::GFFWriter - Take results from non-LTR search and make an annotated GFF

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

has fastadir => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has outdir   => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );

sub write_gff {
    my $self = shift;
    my $outdir = $self->outdir;

    my @all_clade = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1', 'R2', 'CRE');
    my $seq_dir   = File::Spec->catdir($outdir, 'info', 'full');
    my @cladedirs = map { File::Spec->catdir($seq_dir, $_) } @all_clade;
    my @nonltrs;

    for my $clade (@cladedirs) {
        my $name = basename($clade);
	my $filename = $name.".dna";
        if (-d $clade) {
	    find( sub { push @nonltrs, $File::Find::name if -f and /$filename/ }, $clade ); 	    
        }
    }

    $self->_fasta_to_gff(\@nonltrs);
}

sub _fasta_to_gff {
    my $self = shift;
    my $outdir = $self->outdir;
    my ($seqs) = @_;

    my $name = basename($outdir);
    my $outfile = File::Spec->catfile($outdir, $name."_tephra_nonltr.gff3");
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    my $lens = $self->_get_seq_region;

    my %regions;
    for my $file (@$seqs) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	open my $in, '<', $file or die "\nERROR: Could not open file: $file\n";
	while (my $line = <$in>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		my ($id, $start, $end) = ($line =~ /^\>(\w+\.?\w+)_(\d+)-(\d+)/);
		if (defined $id) {
		    $regions{$name}{$id}{$start} = join "||", $start, $end;
		}
	    }
	}
	close $in;
    }

    my @ids = map { nsort keys %{$regions{$_}} } keys %regions;
    say $out '##gff-version 3';

    for my $seqid (@ids) {
	my ($name, $path, $suffix) = fileparse($seqid, qr/\.[^.]*/);
	if (exists $lens->{$name}) {
	    say $out join q{ }, '##sequence-region', $name, '1', $lens->{$name};
	}
	else {
	    say "\nERROR: Could not find $name in map.\n";
	}
    }

    ##TODO: How to get the strand correct?
    my $ct = 0;
    for my $clade (keys %regions) {
	for my $seqid (nsort keys %{$regions{$clade}}) {
	    my ($name, $path, $suffix) = fileparse($seqid, qr/\.[^.]*/);
	    for my $start (sort { $a <=> $b } keys %{$regions{$clade}{$seqid}}) {
	        my ($start, $end) = split /\|\|/, $regions{$clade}{$seqid}{$start};	
	        say $out join "\t", $name, 'Tephra', 'non_LTR_retrotransposon', $start, $end, '-', '?', '-', 
	            "ID=non_LTR_retrotransposon$ct;Name=$clade;Ontology_term=SO:0000189"; 
		$ct++;
	    }
	} 
    }
}

sub _get_seq_region {
    my $self = shift;
    my $fasdir = $self->fastadir;

    my (@seqs, %lens);
    find( sub { push @seqs, $File::Find::name if -f and /\.fa.*$/ }, $fasdir );

    for my $seq (@seqs) {
	my $seqio = Bio::SeqIO->new(-file => $seq, -format => 'fasta');
	while (my $seqobj = $seqio->next_seq) {
	    my $id  = $seqobj->id;
	    my $len = $seqobj->length;
	    $lens{$id} = $len;
	}
    }

    return \%lens;
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

    perldoc Tephra::NonLTR::GFFWriter


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
