package Tephra;

use Moose;
use App::Cmd::Setup -app;
use namespace::autoclean;

=head1 NAME

Tephra - Transposable element palenontology

=head1 VERSION

Version 0.02.5

=cut

our $VERSION = '0.02.5';
#$VERSION = eval $VERSION;

=head1 SYNOPSIS

    tephra <command> [-h] [-m]

=cut

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
