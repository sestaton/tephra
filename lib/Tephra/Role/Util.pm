package Tephra::Role::Util;

use 5.010;
use Moose::Role;
use IPC::System::Simple qw(system);
use Try::Tiny;
use namespace::autoclean;

sub run_cmd {
    my $self = shift;
    my ($cmd) = @_;

    try {
	system([0..5], $cmd);
    }
    catch {
	die "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
}

1;
