package Tephra::Role::Run::Blast;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use Log::Any            qw($log);
use IPC::System::Simple qw(system capture);
use Cwd                 qw(abs_path);
use Try::Tiny;
use File::Spec;
use File::Find;
use File::Basename;
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::Role::Run::Blast - Helper role for running NCBI BLAST

=head1 VERSION

Version 0.09.4

=cut

our $VERSION = '0.09.4';
$VERSION = eval $VERSION;

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
    my ($blastbin) = @{$config}{qw(blastpath)};
    my $makeblastdb = File::Spec->catfile($blastbin, 'makeblastdb');

    try {
	my @makedbout = capture([0..5],"$makeblastdb -in $db_fas -dbtype nucl -title $db -out $db_path 2>&1 > /dev/null");
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
    #my ($dbname, $dbpath, $dbsuffix) = fileparse($db, qr/\.[^.]*/);
    #my ($qname, $qpath, $qsuffix) = fileparse($query, qr/\.[^.]*/);
    #my $blast_report = File::Spec->catfile( abs_path($qpath), $qname."_$dbname".'.bln' );

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($blastbin) = @{$config}{qw(blastpath)};
    my $blastn = File::Spec->catfile($blastbin, 'blastn');
    
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
	my @makedbout = capture([0..5], $cmd);
    }
    catch {
	say STDERR "Unable to run blast. Here is the exception: $_.";
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
