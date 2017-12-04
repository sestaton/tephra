package Tephra::Command::findfragments;
# ABSTRACT: Search a masked genome with a repeat database to find fragmented elements.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use File::Path    qw(make_path remove_tree);
use Tephra -command;
use Tephra::Genome::FragmentSearch;

sub opt_spec {
    return (    
	[ "genome|g=s",     "The masked genome sequences in FASTA format used to search for transposon fragments " ],
	[ "repeatdb|d=s",   "The file of repeat sequences in FASTA format used to query the genome "               ], 
	[ "outfile|o=s",    "The GFF3 file of non-overlapping transposon fragments "                               ],         
	[ "percentid|p=i",  "The percent identity cutoff for BLAST hits to the repeat database (Default: 80) "     ],
	[ "hitlen|l=i",     "The minimum length BLAST hits to the repeat database (Default: 100) "                 ],
	[ "threads|t=i",    "The number of threads to use for parallel BLAST searches (Default: 1) "               ],
    );
}

sub validate_args {
    my ($self, $opt, $args) = @_;

    my $command = __FILE__;
    if ($opt->{man}) {
        system('perldoc', $command) == 0 or die $!;
        exit(0);
    }
    elsif ($opt->{help}) {
        $self->help and exit(0);
    }
    elsif (!$opt->{genome} || !$opt->{repeatdb} || !$opt->{outfile}) {
	say STDERR "\nERROR: Required arguments not given.\n";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _find_transposon_fragments($opt);
}

sub _find_transposon_fragments {
    my ($opt) = @_;

    my $genome    = $opt->{genome};
    my $repeatdb  = $opt->{repeatdb};
    my $outfile   = $opt->{outfile};
    my $hpid      = $opt->{percentid} // 80;
    my $hlen      = $opt->{hitlen} // 100;
    my $threads   = $opt->{threads} // 1;
    
    my $ff_obj = Tephra::Genome::FragmentSearch->new( 
	genome    => $genome, 
	repeatdb  => $repeatdb, 
	outfile   => $outfile, 
	matchlen  => $hlen, 
	percentid => $hpid, 
	threads   => $threads
    );
    $ff_obj->find_transposon_fragments;
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra findfragments [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The masked genome sequences in FASTA format used to search for transposon fragments.
    -d|repeatdb   :   The file of repeat sequences in FASTA format used to query the genome.
    -o|outfile    :   The GFF3 file of non-overlapping transposon fragments.
    
Options:
    -t|threads    :   The number of threads to use for parallel BLAST searches (Default: 1).
    -p|percentid  :   The percent identity cutoff for BLAST alignments to the genome (Default: 80).
    -l|hitlen     :   The minimum length BLAST alignments to the genome (Default: 100).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findfragments -  Search a masked genome with a repeat database to find fragmented elements

=head1 SYNOPSIS    

 tephra findfragments -f custom_repeats.fas -d repeatdb.fas -o ref_classified.fas

=head1 DESCRIPTION

 This subcommand takes a database of full-length transposons and searches a masked genome to identify
 truncated or fragmented elements. The input database can be from Tephra or any other source. If the 'all' command
 is run then the results of this command will added to the final GFF3.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The masked genome sequences in FASTA format used to search for transposon fragments.

=item -d, --repeatdb

 The file of repeat sequences in FASTA format used to query the genome.

=item -o, --outfile

 The GFF3 file of non-overlapping transposon fragments.

=back

=head1 OPTIONS

=over 2

=item -t, --threads

 The number of threads to use for parallel BLAST searches (Default: 1).

=item -p, --percentid

 The percent identity cutoff for BLAST hits to the repeat database (Default: 80).

=item -l, --hitlen

 The minimum length BLAST hits to the repeat database (Default: 100).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
