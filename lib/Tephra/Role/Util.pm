package Tephra::Role::Util;

use 5.014;
use Moose::Role;
use Bio::DB::HTS::Faidx;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Util - Helper methods for running programs

=head1 VERSION

Version 0.09.5

=cut

our $VERSION = '0.09.5';
$VERSION = eval $VERSION;

sub index_ref {
    my $self = shift;
    my ($fasta) = @_;

    my $index = Bio::DB::HTS::Faidx->new($fasta);
    return $index;
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
