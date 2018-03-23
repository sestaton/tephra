package Tephra::Role::Util;

use 5.014;
use Moose::Role;
use File::Basename qw(fileparse);
use File::Temp     qw(tempfile);
#use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Util - Helper methods for running programs

=head1 VERSION

Version 0.09.9

=cut

our $VERSION = '0.09.9';
$VERSION = eval $VERSION;

sub index_ref {
    my $self = shift;
    my ($fasta) = @_;

    # This method is to test whether some IDs cause indexing issues. Though,
    # the hard coded method in place modifies the FASTA names, so this approach should
    # be avoided in practice.
    #my $genome = $self->_adjust_identifiers($fasta);
    my $index = Bio::DB::HTS::Faidx->new($fasta);

    return $index;
}

sub get_full_seq {
    my $self = shift;
    my ($index, $chromosome, $start, $end) = @_;
    my $location = "$chromosome:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);

    return ($seq, $length);
}

sub _adjust_identifiers {
    my $self = shift;
    my ($fasta) = @_;

    my ($name, $path, $suffix) = fileparse($fasta, qr/\.[^.]*/);
    my $tmpiname  = $name.'_XXXX';
    my ($tmp_fh, $tmp_filename) = tempfile( TEMPLATE => $tmpiname, DIR => $path, SUFFIX => $suffix, UNLINK => 0 );

    my $kseq = Bio::DB::HTS::Kseq->new($fasta);
    my $iter = $kseq->iterator;

    while (my $seqo = $iter->next_seq) {
	my $id = $seqo->name;
	my $seq = $seqo->seq;
	$id =~ s/-/_/g;
	say $tmp_fh join "\n", ">$id", $seq;
    }
    close $tmp_fh;

    return $tmp_filename;
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

    perldoc Tephra::Role::Util


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
