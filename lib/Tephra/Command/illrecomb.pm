package Tephra::Command::illrecomb;
# ABSTRACT: Characterize the distribution of illigetimate recombination in a genome.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::Genome::IllRecombination;

sub opt_spec {
    return (
	[ "indir|i=s",       "The directory of LTR families in FASTA format. " ],
	[ "outfile|o=s",     "The filename to write the extracted sequences to " ],
	[ "statsfile|s=s",   "The file to write the alignment stats for all gaps to " ],
	[ "recombstats|r=s", "The file to write the alignment stats for illegetimate recombination sites " ],
	[ "repeatpid|p=i",   "The percent identity threshold for retaining repeats that flank gaps. " ],
	[ "threads|t=i",     "The number of threads to use for alignments (Default: 1) " ],
	[ "clean|c",         "Clean up the intermediate alignment files (Default: yes) " ],
	);
}

sub validate_args {
    my ($self, $opt, $args) = @_;

    my $command = __FILE__;
    if ($self->app->global_options->{man}) {
	system([0..5], "perldoc $command");
    }
    elsif ($self->app->global_options->{help}) {
	$self->help;
    }
    elsif (!$opt->{indir} || !$opt->{outfile} || !$opt->{statsfile}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $some = _calculate_ill_recomb($opt);
}

sub _calculate_ill_recomb {
    my ($opt) = @_;

    my $dir          = $opt->{indir};
    my $outfile      = $opt->{outfile};
    my $allstatsfile = $opt->{statsfile};
    my $illrecstats  = $opt->{recombstats};
    my $threads      = defined $opt->{threads} ? $opt->{threads} : 1;
    my $clean        = defined $opt->{clean} ? 1 : 0;

    my $ill_obj = Tephra::Genome::IllRecombination->new(
	dir             => $dir,
	outfile         => $outfile,
	allstatsfile    => $allstatsfile,
	illrecstatsfile => $illrecstats,
	threads         => $threads,
	clean           => $clean,
    );
    
    $ill_obj->find_illigetimate_recombination;

}

sub help {
    print STDERR<<END

USAGE: tephra illrecomb [-h] [-m]
    -m --man         :   Get the manual entry for a command.
    -h --help        :   Print the command usage.

Required:
    -i|indir         :    The directory of LTR families in FASTA format.
    -o|outfile       :    File name to write the extracted sequences to.

Options:
    -t|threads       :    The number of threads to use for alignments (Default: 1).    
    -r|recombstats   :    The file to write the alignment stats for illegetimate recombination sites.    
    -s|statsfile     :    A file to write alignment stats to.
    -p|repeat_pid    :    The percent identity threshold for retaining repeats that flank gaps. 
                          (Default: 10, which is a 2 bp match).
    -c|clean         :     Clean up the intermediate alignment files (Default: yes).
			  
END
}

=pod

=head1 NAME
                                                                       
 tephra illrecomb

=head1 SYNOPSIS    

 tephra illrecomb

=head1 DESCRIPTION

 This subcommand takes a 

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -d, --repeatdb

 The file of repeat sequences in FASTA format to use for classification.

=item -f, --gff

 The GFF3 file of LTR-RTs in <--genome> as output by the 'tephra findltrs' command.

=item -o, --outdir

 The output directory for placing categorized elements.

=back

=head1 OPTIONS

=over 2

=item -t, --threads

 The number of threads to use for clustering coding domains (Default: 1).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
    
1;
