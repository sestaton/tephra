package Tephra::Role::Run::GT;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use IO::File;
use POSIX;
use Env                 qw(@PATH);
use Cwd                 qw(getcwd abs_path);
use Log::Any            qw($log);
use IPC::System::Simple qw(system EXIT_ANY);
use File::Path          qw(make_path);
use File::Temp;
use File::Spec;
use File::Find;
use File::Basename;
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::Role::Run::GT - Helper role for running genometools

=head1 VERSION

Version 0.09.1

=cut

our $VERSION = '0.09.1';
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
    my ($args, $logfile) = @_;
    my $threads = $self->threads;
    my $log     = $self->get_tephra_logger($logfile);

    my $gt = $self->get_gt_exec;
    unshift @$args, 'suffixerator';
    unshift @$args, "$gt -j $threads";

    my $err = $self->_make_gt_errorlog('suffixerator');
    my $cmd = join qq{ }, @$args;
    $cmd .= " 2> $err";
    say STDERR "DEBUG: $cmd" if $self->debug;
    my $EXIT = system(EXIT_ANY, $cmd);
    my $errors = $self->_errorlog_to_string($err);

    if ($EXIT) { # non-zero
	$log->error("'gt suffixerator' failed with exit value: $EXIT. Here is the output: $errors\n");
	exit(1);
    }

    unlink $err unless -s $err;
    return $EXIT == 0 ? 1 : 0;
}

sub run_ltrharvest {
    my $self = shift;
    my ($args, $logfile) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;
    my $log     = $self->get_tephra_logger($logfile);

    my $gt = $self->get_gt_exec;
    my @ltrh_args;
    push @ltrh_args, "$gt -j $threads ltrharvest";

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrh_args, $opt.q{ }.$args->{$opt};
	}
    } 

    my $err = $self->_make_gt_errorlog('ltrharvest');
    my $cmd = join qq{ }, @ltrh_args;
    $cmd .= " 2>&1 >$err";
    say STDERR "DEBUG: $cmd" if $self->debug; # and exit;
    my $EXIT = system(EXIT_ANY, $cmd);
    my $errors = $self->_errorlog_to_string($err);

    if ($EXIT) { # non-zero  
	$log->error("LTRharvest failed with exit value: $EXIT. Here is the output: $errors\n");
	exit(1);
    }

    unlink $err unless -s $err;
    return $EXIT == 0 ? 1 : 0;
}

sub run_ltrdigest {
    my $self = shift;
    my ($args, $gff, $logfile) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;
    my $log     = $self->get_tephra_logger($logfile);

    # see: https://github.com/genometools/genometools/issues/662
    # there is something wrong with setting the TMPDIR with ltrdirgest
    #my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    #my $time = POSIX::strftime("%d-%m-%Y_%H-%M-%S", localtime);
    #my $model_dir = File::Spec->catdir($path, "tephra_hmm_modeldir_$time");
    #make_path( $model_dir, {verbose => 0, mode => 0771,} );
    #$TMDIR = $model_dir; 

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer3bin) = @{$config}{qw(hmmer3bin)};
    push @PATH, $hmmer3bin;
  
    my $gt = $self->get_gt_exec;
    my @ltrd_args;
    push @ltrd_args, "$gt -j $threads ltrdigest";

    my $hmmdb = $args->{'-hmms'};
    push @ltrd_args, '-hmms'.q{ }.$hmmdb; 
    delete $args->{'-hmms'};

    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @ltrd_args, $opt.q{ }.$args->{$opt};
	}
    }
    
    my $err = $self->_make_gt_errorlog('ltrdigest');
    my $cmd = join qq{ }, @ltrd_args, $gff;
    
    $cmd .= " 2>&1 >$err";
    say STDERR "DEBUG: $cmd" if $self->debug;
    my $EXIT = system(EXIT_ANY, $cmd);
    my $errors = $self->_errorlog_to_string($err);
    undef @PATH;
    #undef $TMPDIR;

    if ($EXIT) { # non-zero 
	$log->error("LTRdigest failed with exit value: $EXIT. Here is the output: $errors\n");
	exit(1);
    }

    unlink $err unless -s $err;
    return $EXIT == 0 ? 1 : 0;
}

sub run_tirvish {
    my $self = shift;
    my ($args, $gff, $logfile) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;
    my $gt      = $self->get_gt_exec;
    my $log     = $self->get_tephra_logger($logfile);

    #my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    #my $time = POSIX::strftime("%d-%m-%Y_%H:%M:%S", localtime);
    #my $model_dir = File::Spec->catdir($path, "tephra_hmm_modeldir_$time");
    #make_path( $model_dir, {verbose => 0, mode => 0771,} );
    #$ENV{TMPDIR} = $model_dir;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer3bin) = @{$config}{qw(hmmer3bin)};
    #$ENV{PATH} = join ':', $ENV{PATH}, $hmmer3bin;
    push @PATH, $hmmer3bin;

    my @tirv_args;
    push @tirv_args, "$gt -j $threads tirvish";
    for my $opt (keys %$args) {
	if (defined $args->{$opt}) {
	    push @tirv_args, $opt.q{ }.$args->{$opt};
	}
    }
    
    my $err = $self->_make_gt_errorlog('tirvish');
    my $cmd = join qq{ }, @tirv_args, ">", $gff;
     $cmd .= " 2> $err";
    say STDERR "DEBUG: $cmd" if $self->debug;
    my $EXIT = system(EXIT_ANY, $cmd);
    my $errors = $self->_errorlog_to_string($err);
    undef @PATH;

    if ($EXIT) { # non-zero 
	$log->error("'gt tirvish' failed exit value: $EXIT. Here is the output: $errors\n");
	exit(1);
    }

    unlink $err unless -s $err;
    return $EXIT == 0 ? 1 : 0;
}

sub sort_gff {
    my $self = shift;
    my ($gff, $logfile) = @_;
    my $debug   = $self->debug;
    my $threads = $self->threads;
    my $log     = $self->get_tephra_logger($logfile);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $gff_sort = File::Spec->catfile( abs_path($path), $name."_sort$suffix" );
    my $gt = $self->get_gt_exec;

    my $err = $self->_make_gt_errorlog('gt-sort');
    my $cmd = "$gt -j $threads gff3 -sort $gff > $gff_sort";
    $cmd .= " 2> $err";
    say STDERR "DEBUG: $cmd" if $self->debug;
    my $EXIT = system(EXIT_ANY, $cmd);
    my $errors = $self->_errorlog_to_string($err);

    if ($EXIT) { # non-zero 
	$log->error("'gt gff3 -sort' failed exit value: $EXIT. Here is the output: $errors\n");
	exit(1);
    }

    unlink $err unless -s $err;
    return $EXIT == 0 ? $gff_sort : 0;
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

sub _make_gt_errorlog {
    my $self = shift;
    my ($cmd) = @_;
    
    my $dir = getcwd();
    my $tmpiname  = "tephra_$cmd"."_errors_XXXX";
    my $tmp_hmmdbfh = File::Temp->new( TEMPLATE => $tmpiname,
				       DIR      => $dir,
				       SUFFIX   => '.err',
				       UNLINK   => 0);
    my $tmp_hmmdb = $tmp_hmmdbfh->filename;
    
    
    return $tmp_hmmdb;
}

sub _errorlog_to_string {
    my $self = shift;
    my ($log) = @_;

    my $lines = do {
	local $/ = undef;
	open my $in, '<', $log or die "\nERROR: Could not open file: $log\n";
	<$in>;
    };
    unlink $log;

    return $lines;
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
        #$log->error("Unable to find 'gt' executable. Make sure genometools is installed. Exiting.");
	say STDERR "\nERROR: Unable to find 'gt' executable. Make sure genometools is installed. Exiting.\n";
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
