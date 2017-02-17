package Tephra::Role::Run::GT;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use Cwd                 qw(getcwd abs_path);
use Log::Any            qw($log);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use File::Spec;
use File::Find;
use File::Basename;
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::GT - Helper role for running genometools

=head1 VERSION

Version 0.06.1

=cut

our $VERSION = '0.06.1';
$VERSION = eval $VERSION;

has gt_exec => (
    is        => 'rw',
    isa       => 'Path::Class::File',
    required  => 0,
    coerce    => 1, 
    reader    => 'get_gt_exec',
    writer    => 'set_gt_exec',
    predicate => 'has_gt_exec',
    builder   => '_build_gt_exec',
);

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has hmmdb => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 0,
    coerce   => 1,
);

has trnadb => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 0,
    coerce   => 1,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

has debug => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

sub create_index {
    my $self = shift;
    my ($args) = @_;

    my $gt = $self->get_gt_exec;
    unshift @$args, 'suffixerator';
    unshift @$args, $gt;
    my $cmd = join qq{ }, @$args;

    say STDERR "DEBUG: $cmd" if $self->debug;
    my ($stdout, $stderr, $exit);
    try {
	system($cmd);
	#my @out = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("Unable to make suffixerator index. Here is the exception: $_\nExiting.");
	exit(1);
    };
}

sub run_ltrharvest {
    my $self = shift;
    my ($args) = @_;
    my $debug = $self->debug;

    my $gt = $self->get_gt_exec;
    my @ltrh_args;
    push @ltrh_args, "$gt ltrharvest";

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrh_args, $opt.q{ }.$args->{$opt};
	}
    } 

    my $cmd = join qq{ }, @ltrh_args;
    say STDERR "DEBUG: $cmd" if $self->debug; # and exit;
    try {
	#system($cmd);
	my @out = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("LTRharvest failed. Here is the exception: $_\nExiting.");
	exit(1);
    };
}

sub run_ltrdigest {
    my $self = shift;
    my ($args, $gff) = @_;
    my $debug = $self->debug;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer3bin) = @{$config}{qw(hmmer3bin)};
    $ENV{PATH} = join ':', $ENV{PATH}, $hmmer3bin;

    my $gt = $self->get_gt_exec;
    my @ltrd_args;
    push @ltrd_args, "$gt ltrdigest";
    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrd_args, $opt.q{ }.$args->{$opt};
	}
    }
    
    my $cmd = join qq{ }, @ltrd_args, $gff;
    say STDERR "DEBUG: $cmd" if $self->debug;
    try {
	#system($cmd);
	my @out = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("LTRdigest failed. Here is the exception: $_\nExiting.");
	exit(1);
    };
}

sub run_tirvish {
    my $self = shift;
    my ($args, $gff) = @_;
    my $debug = $self->debug;
    my $gt = $self->get_gt_exec;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer3bin) = @{$config}{qw(hmmer3bin)};
    $ENV{PATH} = join ':', $ENV{PATH}, $hmmer3bin;

    my @tirv_args;
    push @tirv_args, "$gt tirvish";
    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @tirv_args, $opt.q{ }.$args->{$opt};
	}
    }

    my $cmd = join qq{ }, @tirv_args, ">", $gff;
    say STDERR "DEBUG: $cmd" if $self->debug;
    try {
	my @out = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("'gt tirvish' failed. Here is the exception: $_\nExiting.");
	exit(1);
    };
}

sub sort_gff {
    my $self = shift;
    my ($gff) = @_;
    my $debug = $self->debug;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $gff_sort = File::Spec->catfile( abs_path($path), $name."_sort".$suffix );
    my $gt = $self->get_gt_exec;

    my $cmd = "$gt gff3 -sort $gff > $gff_sort";
    say STDERR "DEBUG: $cmd" if $self->debug;
    try {
	my @out = capture { system([0..5], $cmd) };
    }
    catch {
	$log->error("'gt gff3 -sort' failed. Here is the exception: $_\nExiting.");
	exit(1);
    };

    return $gff_sort;
}

sub clean_indexes {
    my $self = shift;
    my ($dir) = @_;

    my @files;
    find( sub { push @files, $File::Find::name
		    if /\.llv|\.md5|\.prf|\.tis|\.suf|\.lcp|\.ssp|\.sds|\.des|\.dna|\.esq|\.prj|\.ois/
	  }, $dir);
    unlink @files;
}

sub clean_index_files {
    my $self = shift;
    my ($index) = @_;
    
    my ($name, $path, $suffix) = fileparse($index, qr/\.[^.]*/);

    my $pat;
    for my $suf (qw(.al1 .llv .ssp .bck .ois .sti1 .bwt .prj .suf .des .sds .tis .lcp .skp)) {
        $pat .= "$name$suffix$suf|";
    }
    $pat =~ s/\|$//;

    my @files;
    find( sub { push @files, $File::Find::name if /$pat/ }, $path);

    unlink @files;
    return;
}

sub _build_gt_exec {
    my $self  = shift;
    my $gtexe = $self->get_gt_exec; # set during class initialization
    
    # check if the executable path was set first
    if (defined $gtexe && -e $gtexe && -x $gtexe) {
        return $gtexe;
    }
    elsif (! defined $gtexe) {
	my $config = Tephra::Config::Exe->new->get_config_paths;
	($gtexe) = @{$config}{qw(gt)};
	if (-e $gtexe && -x $gtexe) {
	    return $gtexe;
	}

	my @path = split /:|;/, $ENV{PATH};
	for my $p (@path) {
	    my $gt = File::Spec->catfile($p, 'gt');

	    if (-e $gt && -x $gt) {
		$self->set_gt_exec($gt);
		return $gt;
	    }
	}
    }
    else {
        $log->error("Unable to find 'gt' executable. Make sure genometools is installed. Exiting.");
        exit(1);
    }
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

    perldoc Tephra::Role::Run::GT


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
