package Tephra::Role::Run::GT;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use IO::File;
use Cwd                 qw(getcwd abs_path);
use Log::Any            qw($log);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(capture);
use Try::Tiny;
use File::Spec;
use File::Find;
use File::Basename;
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::GT - Helper role for running genometools

=head1 VERSION

Version 0.08.0

=cut

our $VERSION = '0.08.0';
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

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

sub create_index {
    my $self = shift;
    my ($args) = @_;
    my $threads = $self->threads;

    my $gt = $self->get_gt_exec;
    unshift @$args, 'suffixerator';
    unshift @$args, "$gt -j $threads";
    my $cmd = join qq{ }, @$args;

    say STDERR "DEBUG: $cmd" if $self->debug;
    my ($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    return $exit == 0 ? 1 : 0;

    if ($exit) { # non-zero
	$log->error("'gt suffixerator' failed with exit code: $exit. Here is STDOUT: $stdout\nHere is stderr: $stderr\n");
	exit(1);
    }
}

sub run_ltrharvest {
    my $self = shift;
    my ($args) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;


    my $gt = $self->get_gt_exec;
    my @ltrh_args;
    push @ltrh_args, "$gt -j $threads ltrharvest";
    #push @ltrh_args, $gt;
    #push @ltrh_args, 'ltrharvest';

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrh_args, $opt.q{ }.$args->{$opt};
	}
    } 

    my $cmd = join qq{ }, @ltrh_args;
    say STDERR "DEBUG: $cmd" if $self->debug; # and exit;
    my ($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    return $exit == 0 ? 1 : 0;

    if ($exit) { # non-zero  
	$log->error("LTRharvest failed with exit code: $exit. Here is STDOUT: $stdout\nHere is stderr: $stderr\n");
	exit(1);
    };
}

sub run_ltrdigest {
    my $self = shift;
    my ($args, $gff) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer3bin) = @{$config}{qw(hmmer3bin)};
    $ENV{PATH} = join ':', $ENV{PATH}, $hmmer3bin;

    my $gt = $self->get_gt_exec;
    my @ltrd_args;
    push @ltrd_args, "$gt -j $threads ltrdigest";
    #push @ltrd_args, $gt;
    #push @ltrd_args, 'ltrdigest';

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrd_args, $opt.q{ }.$args->{$opt};
	}
    }
    
    #say STDERR "DEBUG: $gff";
    my $cmd = join qq{ }, @ltrd_args, $gff;
    say STDERR "DEBUG: $cmd" if $self->debug;
    my ($stdout, $stderr, $exit) = capture { system($cmd) };
    #say STDERR "DEBUG: ltrdigest exit status: $exit" if $self->debug;
    return $exit == 0 ? 1 : 0;

    if ($exit) { # non-zero 
	$log->error("LTRdigest failed with exit code: $exit. Here is STDOUT: $stdout\nHere is stderr: $stderr\n");
	exit(1);
    }
}

sub run_tirvish {
    my $self = shift;
    my ($args, $gff) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;
    my $gt      = $self->get_gt_exec;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer3bin) = @{$config}{qw(hmmer3bin)};
    $ENV{PATH} = join ':', $ENV{PATH}, $hmmer3bin;

    my @tirv_args;
    push @tirv_args, "$gt -j $threads tirvish";
    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @tirv_args, $opt.q{ }.$args->{$opt};
	}
    }

    my $cmd = join qq{ }, @tirv_args, ">", $gff;
    say STDERR "DEBUG: $cmd" if $self->debug;
    my ($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    return $exit == 0 ? 1 : 0;

    if ($exit) { # non-zero 
	$log->error("'gt tirvish' failed exit code: $exit. Here is STDOUT: $stdout\nHere is stderr: $stderr\n");
	exit(1);
    }
}

sub sort_gff {
    my $self = shift;
    my ($gff) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $gff_sort = File::Spec->catfile( abs_path($path), $name."_sort".$suffix );
    my $gt = $self->get_gt_exec;

    my $cmd = "$gt -j $threads gff3 -sort $gff > $gff_sort";
    say STDERR "DEBUG: $cmd" if $self->debug;
    my ($stdout, $stderr, $exit) = capture { system([0..5], $cmd) };
    return $exit == 0 ? $gff_sort : 0;

    if ($exit) { # non-zero 
	$log->error("'gt gff3 -sort' failed exit code: $exit. Here is STDOUT: $stdout\nHere is stderr: $stderr\n");
	exit(1);
    }
}

sub clean_indexes {
    my $self = shift;
    my ($dir) = @_;

    my @files;
    my $wanted  = sub { push @files, $File::Find::name
			    if -f && /\.llv|\.md5|\.prf|\.tis|\.suf|\.lcp|\.ssp|\.sds|\.des|\.dna|\.esq|\.prj|\.ois/ };
    my $process = sub { grep ! -d, @_ };
    find({ wanted => $wanted, preprocess => $process }, $dir);
    unlink @files;

    return;
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
