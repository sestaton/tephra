package Tephra::Role::GT;

use 5.010;
use Moose::Role;
use MooseX::Types::Path::Class;
use IPC::System::Simple qw(system);
use Try::Tiny;
use File::Spec;
use File::Find;
use File::Basename;
use Log::Any        qw($log);
use Cwd;
use namespace::autoclean;

use Data::Dump;

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

sub create_index {
    my $self = shift;
    my ($args) = @_;

    my $gt = $self->get_gt_exec;
    ##TODO
    #my @ltrh_args;
    unshift @$args, 'suffixerator';
    unshift @$args, $gt;
    my $cmd = join qq{ }, @$args;
    #for my $opt (keys %$args) {
	#push @ltrh_args, $opt.q{ }.$args->{$opt};
    #}
    #my $index = File::Spec->catfile($path, $genome.".index");
    #push @ltrh_args, "-indexname $index";
    #my $index_cmd = "$gt suffixerator ";
    #$index_cmd .= "-db $genome -indexname $index -tis -suf -lcp -ssp -sds -des -dna";
    #say STDERR join qq{ }, $cmd;
    try {
	system([0..5], $cmd);
    }
    catch {
	$log->error("Unable to make suffixerator index. Here is the exception: $_\nExiting.");
	exit(1);
    };

    #return $index;
}

sub run_ltrharvest {
    my $self = shift;
    my ($args) = @_;

    #dd $args and exit;
    my $gt = $self->get_gt_exec;
    my @ltrh_args;
    push @ltrh_args, "$gt ltrharvest";
    #push @ltrh_args, 'ltrharves
    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrh_args, $opt.q{ }.$args->{$opt};
	}
    } 
    #my $ltrh_cmd = "$gt ltrharvest $args 2>&1 > /dev/null";
    #say STDERR $ltrh_cmd;
    #unshift @$args, 'ltrharvest';
    #unshift @$args, $gt;
    my $cmd = join qq{ }, @ltrh_args;
    #say STDERR $cmd and exit;
    try {
	system([0..5], $cmd);
    }
    catch {
	$log->error("LTRharvest failed. Here is the exception: $_\nExiting.");
	exit(1);
    };

}

sub run_ltrdigest {
    my $self = shift;
    my ($args, $gff) = @_;

    my $gt = $self->get_gt_exec;
    #my $ltrd_cmd = "$gt ltrdigest $args";
    my @ltrd_args;
    push @ltrd_args, "$gt ltrdigest";
    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrd_args, $opt.q{ }.$args->{$opt};
	}
    }
    
    #unshift @$args, 'ltrdigest';
    #unshift @$args, $gt;

    my $cmd = join qq{ }, @ltrd_args, $gff;
    #say STDERR $cmd and exit;
    try {
	system([0..5], $cmd);
    }
    catch {
	$log->error("LTRdigest failed. Here is the exception: $_\nExiting.");
	exit(1);
    };

}

sub run_tirvish {
    my $self = shift;
    my ($args) = @_;

    my $hmms  = $args->{hmms};
    my $index = $args->{index};
    my $gff   = $args->{gff};
    my $gt    = $self->get_gt_exec;
    my $cmd = "$gt tirvish -index $index -hmms $hmms > $gff";
    try {
	system([0..5], $cmd);
    }
    catch {
	$log->error("'gt tirvish' failed. Here is the exception: $_\nExiting.");
	exit(1);
    };

    #return $gff;
}

sub sort_gff {
    my $self = shift;
    my ($gff) = @_;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $gff_sort = File::Spec->catfile($path, $name."_sort".$suffix);
    my $gt = $self->get_gt_exec;

    my $sort_cmd = "$gt gff3 -sort $gff > $gff_sort";
    try {
	system([0..5], $sort_cmd);
    }
    catch {
	$log->error("'gt gff3 -sort' failed. Here is the exception: $_\nExiting.");
	exit(1);
    };

    return $gff_sort;
}

sub clean_index {
    my $self = shift;

    my $dir = getcwd();
    my @files;
    find( sub { push @files, $File::Find::name
		    if /\.llv|\.md5|\.prf|\.tis|\.suf|\.lcp|\.ssp|\.sds|\.des|\.dna|\.esq|\.prj|\.ois/
	  }, $dir);
    unlink @files;
}

sub _build_gt_exec {
    my $self  = shift;
    my $gtexe = $self->get_gt_exec; # set during class initialization
    #my $gtexe = '/home';
    
    # check if the executable path was set first
    if (defined $gtexe && -e $gtexe && -x $gtexe) {
        return $gtexe;
    }
    elsif (! defined $gtexe) {
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
        #$log->error(
	say STDERR "Unable to find 'gt' executable. Make sure genometools is installed. Exiting."; #);
        exit(1);
    }
}

1;
