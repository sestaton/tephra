package Tephra::Role::Run::HelitronScanner;

use 5.010;
use Moose::Role;
use MooseX::Types::Path::Class;
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use File::Spec;
use File::Find;
use File::Basename;
use Log::Any        qw($log);
use Cwd;
use namespace::autoclean;

#use Data::Dump;
#use Data::Printer;

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has helitronscanner_dir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 1,
    coerce   => 1,
);

has outfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

sub run_hscan_headtail {
    my $self = shift;
    my ($args, $jar, $subcmd) = @_;

    #dd $args and exit;
    my @scan_args;
    push @scan_args, "java -jar $jar $subcmd";

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @scan_args, $opt.q{ }.$args->{$opt};
	}
    } 
    push @scan_args, "--rc";
    
    my $cmd = join qq{ }, @scan_args;
    say STDERR $cmd;# and exit;

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
    say STDERR $cmd; # and exit;
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
    say STDERR $cmd; # and exit;
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

1;
