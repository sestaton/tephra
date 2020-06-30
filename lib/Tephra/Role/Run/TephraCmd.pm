package Tephra::Role::Run::TephraCmd;

use 5.014;
use Moose::Role;
use IPC::System::Simple qw(system capture);
use Try::Tiny;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::TephraCmd - Helper role for running Tephra subcommands

=head1 VERSION

Version 0.13.0

=cut

our $VERSION = '0.13.0';
$VERSION = eval $VERSION;

sub run_tephra_cmd {
    my $self = shift;
    my ($subcmd, $opts, $debug) = @_;

    say STDERR join q{ }, 'DEBUG: ', 'tephra', $subcmd, @$opts 
        if $debug;

    try {
        my @run_out = system([0..5], 'tephra', $subcmd, @$opts);
    }
    catch {
        say STDERR "Unable to run 'tephra $subcmd'. Here is the exception: $_.";
        exit(1);
    };
}

sub capture_tephra_cmd {
    my $self = shift;
    my ($subcmd, $opts, $debug) = @_;

    say STDERR join q{ }, 'DEBUG: ', 'tephra', $subcmd, @$opts
        if $debug;

    try {
        my @run_out = capture([0..5], 'tephra', $subcmd, @$opts);
    }
    catch {
        say STDERR "Unable to run 'tephra $subcmd'. Here is the exception: $_.";
        exit(1);
    };
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

    perldoc Tephra::Role::Run::TephraCmd


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
