package Tephra::Role::Run::Any;

use 5.014;
use Moose::Role;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::Any - Helper role for running shell commands

=head1 VERSION

Version 0.12.6

=cut

our $VERSION = '0.12.6';
$VERSION = eval $VERSION;

sub run_cmd {
    my $self = shift;
    my ($cmd) = @_;

    try {
        system([0..5], $cmd);
    }
    catch {
        die "\n[ERROR]: $cmd failed. Here is the exception: $_\n";
    };
}

sub capture_cmd {
    my $self = shift;
    my ($cmd) = @_;

    try {
        my @out = capture { system([0..5], $cmd) };
    }
    catch {
        my $err = $_;
        if ($err =~ /SEGV/) {
            return 'failed';
        }
        die "\n[ERROR]: $cmd failed. Here is the exception: $_\n";
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

    perldoc Tephra::Role::Run::Any


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
