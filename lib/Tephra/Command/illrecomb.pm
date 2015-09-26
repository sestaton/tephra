package Tephra::Command::illrecomb;
# ABSTRACT: Characterize the distribution of illigetimate recombination in a genome.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::Genome::IllRecombination;

sub opt_spec {
    return (
	[ "infile|i=s",    "An alignment file in FASTA format. " ],
	[ "outfile|o=s",   "The filename to write the extracted sequences to " ],
	[ "statsfile=s",   "The file to write the alignment stats to " ],
	[ "repeatpid|p=i", "The percent identity threshold for retaining repeats that flank gaps. " ],
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
    elsif (!$opt->{genome} || !$opt->{gff}) {
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


}

sub help {
    print STDERR<<END

USAGE: tephra illrecomb [-h] [-m]
    -m --man         :   Get the manual entry for a command.
    -h --help        :   Print the command usage.

Required:
    -i|infile        :    An alignment file in FASTA format.
    -o|outfile       :    File name to write the extracted sequences to.

Options:
    -s|statsfile     :    A file to write alignment stats to.
    -p|repeat_pid    :    The percent identity threshold for retaining repeats that flank gaps. 
                          (Default: 10, which is a 2 bp match).

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
