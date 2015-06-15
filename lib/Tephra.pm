package Tephra;

use Moo;
use YAML::Tiny;
use Log::Any qw($log);
use namespace::clean;

with 'Tephra::Role::Config';

=head1 NAME

Tephra - Transposable element palenontology

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
#$VERSION = eval $VERSION;

=head1 SYNOPSIS

    tephra --config tephra_config.yml

=cut

has 'config' => (
    is            => 'ro',
    isa           => 'Str',
    required      => 0,
    documentation => qq{The Tephra configuration file},
);

sub get_configuration {
    my $self = shift;
    my $configfile   = YAML::Tiny->read( $self->config );
    my $valid_config = $self->parse_configuration( $configfile );
    return $valid_config;
}

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/Tephra/issues>. I will be notified,
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

1;
