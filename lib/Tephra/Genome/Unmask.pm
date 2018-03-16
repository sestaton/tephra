package Tephra::Genome::Unmask;

use 5.014;
use Moose;
use File::Spec;
use File::Basename;
use Cwd         qw(abs_path);
use File::Copy  qw(move);
use Tephra::Annotation::Util;
use namespace::autoclean;
#use Data::Dump::Color qw(dump dd);

with 'Tephra::Role::Util';

=head1 NAME

Tephra::Genome::Unmask - Unmask a masked sequence using a reference genome 

=head1 VERSION

Version 0.09.8

=cut

our $VERSION = '0.09.8';
$VERSION = eval $VERSION;

has genome => (
      is       => 'ro',
      isa      => 'Maybe[Str]',
      required => 1,
      coerce   => 0,
);

has repeatdb => (
      is       => 'ro',
      isa      => 'Maybe[Str]',
      required => 0,
      coerce   => 0,
);

has outfile => (
    is       => 'ro',
    isa      => 'Maybe[Str]',
    required => 0,
    coerce   => 0,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

sub unmask_repeatdb {
    my $self = shift;
    my $genome   = $self->genome;
    my $repeatdb = $self->repeatdb;
    my $outfile  = $self->outfile // $repeatdb;

    my $index = $self->index_ref($genome);
    my $store = $self->store_seq_coords($repeatdb);

    my ($name, $path, $suffix) = fileparse($repeatdb, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }

    my $tmp_outfile  = File::Spec->catfile( abs_path($path), $name.'_tmp_unmasked.fas' );

    open my $out, '>', $tmp_outfile or die "\nERROR: Could not open file: $tmp_outfile\n";

    for my $id (keys %$store) {
	my $seq = $store->{$id}{seq};
	my ($chr, $start, $end) = split /\|\|/, $store->{$id}{coords};
	my $gseq = $self->get_full_seq($index, $chr, $start, $end);
	$gseq =~ s/.{60}\K/\n/g;

	say $out join "\n", ">".$id, $gseq;
    }
    close $out;

    move $tmp_outfile, $outfile or die "\nERROR: move failed: $!\n";

    return;
}

sub store_seq_coords {
    my $self = shift;
    my ($file) = @_;

    my %hash;
    my $kseq = Bio::DB::HTS::Kseq->new($file);
    my $iter = $kseq->iterator();

    while (my $seqobj = $iter->next_seq) {
        my $id = $seqobj->name;
        #my ($chr, $start, $end) = ($id =~ /_([A-Za-z0-9])_(\d+)_(\d+)$/);
        $id =~ s/_[+-]$//;
        my ($chr, $start, $end) = (split /\_/, $id)[-3..-1];# say join q{ }, $chr, $start, $end'
        my $coords = join "||", $chr, $start, $end;
        my $seq = $seqobj->seq;
        $hash{$id} = { seq => $seq, coords => $coords };
    }

    return \%hash;
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

    perldoc Tephra::Genome::Unmask


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
