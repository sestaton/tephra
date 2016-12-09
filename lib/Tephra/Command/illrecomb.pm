package Tephra::Command::illrecomb;
# ABSTRACT: Characterize the distribution of illegitimate recombination in a genome.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::Genome::IllRecombination;

sub opt_spec {
    return (
	[ "infile|i=s",      "The combined file of LTR families in FASTA format. "                         ],
	[ "outfile|o=s",     "The filename to write the extracted sequences to "                           ],
	[ "statsfile|s=s",   "The file to write the alignment stats for all gaps to "                      ],
	[ "recombstats|r=s", "The file to write the alignment stats for illegetimate recombination sites " ],
	[ "repeatpid|p=i",   "The percent identity threshold for retaining repeats that flank gaps. "      ],
	[ "threads|t=i",     "The number of threads to use for alignments (Default: 1) "                   ],
	[ "clean|c",         "Clean up the intermediate alignment files (Default: yes) "                   ],
	[ "help|h",          "Display the usage menu and exit. "                                           ],
        [ "man|m",           "Display the full manual. "                                                   ],
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
        $self->help;
        exit(0);
    }
    elsif (!$opt->{infile} || !$opt->{outfile} || !$opt->{statsfile}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _calculate_ill_recomb($opt);
}

sub _calculate_ill_recomb {
    my ($opt) = @_;

    my $infile       = $opt->{infile};
    my $outfile      = $opt->{outfile};
    my $allstatsfile = $opt->{statsfile};
    my $illrecstats  = $opt->{recombstats};
    my $threads      = $opt->{threads} // 1;
    my $clean        = $opt->{clean} // 0;

    my $ill_obj = Tephra::Genome::IllRecombination->new(
	infile          => $infile,
	outfile         => $outfile,
	allstatsfile    => $allstatsfile,
	illrecstatsfile => $illrecstats,
	threads         => $threads,
	clean           => $clean,
    );
    
    $ill_obj->find_illegitimate_recombination;
}

sub help {
    print STDERR<<END

USAGE: tephra illrecomb [-h] [-m]
    -m --man         :   Get the manual entry for a command.
    -h --help        :   Print the command usage.

Required:
    -i|infile        :    The combined file of LTR families in FASTA format produced by the
                          'tephra classifyltrs' command.
    -o|outfile       :    File name to write the extracted sequences to.
    -s|statsfile     :    A file to write alignment stats to.

Options:
    -t|threads       :    The number of threads to use for alignments (Default: 1).    
    -r|recombstats   :    The file to write the alignment stats for illegetimate recombination sites.    
    -p|repeat_pid    :    The percent identity threshold for retaining repeats that flank gaps. 
                          (Default: 10, which is a 2 bp match).
    -c|clean         :    Clean up the intermediate alignment files (Default: yes).

END
}

=pod

=head1 NAME
                                                                       
 tephra illrecomb - Calculate the nature and extent of illegitimate recombination in a genome

=head1 SYNOPSIS    

 tephra illrecomb -i tephra_combined_ltrs.fasta -s illrecomb_stats.txt -o illrecomb_seqs.fasta

=head1 DESCRIPTION

 This subcommand calculates the nature of illegitimate recombination (total events, size, sequences involved),
 in a genome by analyzing one TE family at a time. By comparison, it is possible to see if there are
 differences between TE families (my analyses suggest there should not be).

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --infile

 The file of combined LTR families in FASTA format that is produced by the 'tephra classifyltrs' command.

=item -o, --outfile

 File name to write the extracted sequences to.

=item -s, --statsfile

 A file to write alignment stats to.

=head1 OPTIONS

=over 2

=item -t, --threads

 The number of threads to use for alignments (Default: 1).

=item -r, --recombstats

 The file to write the alignment stats for illegetimate recombination sites.

=item -p, --repeat_pid

 The percent identity threshold for retaining repeats that flank gaps.
 The default is 10, which is a 2 bp match.

=item -c, --clean

 Clean up the intermediate alignment files (Default: yes).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
    
1;
