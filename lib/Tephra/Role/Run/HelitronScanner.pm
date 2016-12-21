package Tephra::Role::Run::HelitronScanner;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use Log::Any            qw($log);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use File::Spec;
use File::Find;
use File::Basename;
use Cwd;
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::Role::Run::HelitronScanner - Helper role for running HelitronScanner

=head1 VERSION

Version 0.04.5

=cut

our $VERSION = '0.04.5';
$VERSION = eval $VERSION;

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has helitronscanner => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has gff => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has fasta => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has debug => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

sub run_hscan_headtail {
    my $self = shift;
    my ($args, $jar, $subcmd) = @_;
    
    my @scan_args;
    push @scan_args, "java -jar $jar $subcmd";

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @scan_args, $opt.q{ }.$args->{$opt};
	}
    } 
    push @scan_args, "--rc";
    
    my $cmd = join qq{ }, @scan_args;
    say STDERR $cmd if $self->debug;

    my ($stdout, $stderr, $exit);
    try {
	my @out = capture { system([0..5], $cmd) };
	#($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("HelitronScanner failed. Error: $stderr. Here is the exception: $_\nExiting.");
	exit(1);
    };

}

sub run_hscan_pair {
    my $self = shift;
    my ($args, $jar) = @_;

    my @pair_args;
    push @pair_args, "java -jar $jar pairends";

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @pair_args, $opt.q{ }.$args->{$opt};
	}
    }
    push @pair_args, "--rc";

    my $cmd = join qq{ }, @pair_args;
    say STDERR $cmd if $self->debug;
    my ($stdout, $stderr, $exit);
    try {
	my @out = capture { system([0..5], $cmd) };
	#($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("HelitronScanner failed. Error: $stderr. Here is the exception: $_\nExiting.");
	exit(1);
    };
}

sub run_hscan_draw {
    my $self = shift;
    my ($args, $jar) = @_;

    my @draw_args;
    push @draw_args, "java -jar $jar draw";

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @draw_args, $opt.q{ }.$args->{$opt};
	}
    }

    my $cmd = join qq{ }, @draw_args;
    $cmd .= " --pure --flanking --ext";
    say STDERR $cmd if $self->debug;
    my ($stdout, $stderr, $exit);
    try {
	my @out = capture { system([0..5], $cmd) };
	#($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("HelitronScanner failed. Error: $stderr. Here is the exception: $_\nExiting.");
	exit(1);
    };
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

    perldoc Tephra::Role::Run::HelitronScanner


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
