package Tephra::Role::Run::Blast;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use IPC::System::Simple qw(system capture);
use Cwd                 qw(abs_path);
use Try::Tiny;
use File::Spec;
use File::Basename;
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::Blast - Helper role for running NCBI BLAST

=head1 VERSION

Version 0.13.0

=cut

our $VERSION = '0.13.0';
$VERSION = eval $VERSION;

has infile => (
    is       => 'rw',
    isa      => 'Path::Class::File',
    reader   =>'get_infile',
    writer   =>'set_infile',
    required => 0,
    coerce   => 1,
);

has outfile => (
    is       => 'rw',
    isa      => 'Path::Class::File',
    reader   => 'get_outfile',
    writer   => 'set_outfile',
    required => 0,
    coerce   => 1,
);

has repeatdb => (
    is       => 'rw',
    isa      => 'Path::Class::File',
    reader   => 'get_repeatdb',
    writer   => 'set_repeatdb',
    required => 0,
    coerce   => 1,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

has blast_hit_pid => (
    is      => 'ro',
    isa     => 'Int',
    default => 80,
);

has blast_hit_cov => (
    is      => 'ro',
    isa     => 'Int',
    default => 50,
);

has blast_hit_len => (
    is      => 'ro',
    isa     => 'Int',
    default => 80,
);

sub process_blast_args {
    my $self = shift;
    my $fasta   = $self->get_infile;
    my $rdb     = $self->get_repeatdb;
    my $out     = $self->get_outfile;
    my $threads = $self->get_threads;

    my ($dbname, $dbpath, $dbsuffix) = fileparse($rdb, qr/\.[^.]*/);
    my ($faname, $fapath, $fasuffix) = fileparse($fasta, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($fapath, $faname.'_'.$dbname.'.bln');

    my $blastdb = $self->make_blastdb($rdb);
    my $blast_report = $self->run_blast({ query   => $fasta, 
					  db      => $blastdb, 
					  threads => $threads, 
					  outfile => $outfile,
					  sort    => 'bitscore' });
    my @dbfiles = glob "$blastdb*";
    unlink @dbfiles;

    return $blast_report;
}

sub make_blastdb {
    my $self = shift;
    my ($db_fas) = @_;
    my ($dbname, $dbpath, $dbsuffix) = fileparse($db_fas, qr/\.[^.]*/);

    my $db = $dbname.'_blastdb';
    #my $dir = getcwd();
    my $db_path = Path::Class::File->new( abs_path($dbpath), $db );
    unlink $db_path if -e $db_path;
    #say STDERR "DB: $db_path";

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my $makeblastdb = $config->{makeblastdb};

    try {
	my @makedbout = capture([0..5], "$makeblastdb -in $db_fas -dbtype nucl -title $db -out $db_path 2>&1 > /dev/null");
    }
    catch {
	say STDERR "Unable to make blast database. Here is the exception: $_.";
	say STDERR "Ensure you have removed non-literal characters (i.e., "*" or "-") in your repeat database file.";
	say STDERR "These cause problems with BLAST+. Exiting.";
	exit(1);
    };

    return $db_path;
}

sub run_blast {
    my $self = shift;
    my ($args) = @_;
    my ($query, $db, $threads, $sort, $outfile, $evalue) = @{$args}{qw(query db threads sort outfile evalue)};
    $evalue //= 10;

    # make sure sort can be found under different shells
    my @path = split /:|;/, $ENV{PATH};
    unless (@path) { 
	$ENV{PATH} = "/usr/bin:/bin";
    }

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my $blastn = $config->{blastn};

    my $cmd = "$blastn -query $query -db $db -evalue $evalue -outfmt 6 -num_threads $threads";
    if (defined $sort) {
	if ($sort eq 'bitscore') {
	    $cmd .= " | sort -nrk12,12 >$outfile";
	}
	elsif ($sort eq 'coordinate') {
	    $cmd .= " | sort -nk9,9 >$outfile";
	}
    }
    else {
	$cmd .= " -out $outfile";
    }

    try {
	my @runout = capture([0..5], $cmd);
    }
    catch {
	say STDERR "[ERROR]: Unable to run blast. Here is the exception: $_.";
	exit(1);
    };

    return $outfile;
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

    perldoc Tephra::Role::Run::Blast


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
