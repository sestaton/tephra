package Tephra::Command::tirage;
# ABSTRACT: Calculate the age distribution of TIR transposons.

use 5.014;
use strict;
use warnings;
use Tephra -command;
use Tephra::TIR::TIRStats;

sub opt_spec {
    return (
	[ "genome|g=s",    "The genome sequences in FASTA format used to search for TIRs "              ],
	[ "gff|f=s",       "The GFF3 file of TIRs in <genome> "                                         ],
	[ "outfile|o=s",   "The output file containing the age of each element "                        ],
	[ "subs_rate|r=f", "The nucleotide substitution rate to use (Default: 1e-8) "                   ],
	[ "threads|t=i",   "The number of threads to use for clustering coding domains "                ],
	[ "indir|i=s",     "The input directory of superfamily exemplars "                              ],
	[ "all|a",         "Calculate age of all TIRs in <gff> instead of exemplars in <indir> "        ],
	[ "clean|c",       "Clean up all the intermediate files from PAML and clustalw (Default: yes) " ],
	[ "help|h",        "Display the usage menu and exit. "                                          ],
        [ "man|m",         "Display the full manual. "                                                  ],
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
    elsif (! $opt->{genome} || ! -e $opt->{genome}) {
        say STDERR "\nERROR: The '--genome' file does not appear to exist. Check input.";
        $self->help and exit(0);
    }
    elsif (! $opt->{outfile}) {
	say STDERR "\nERROR: The '--outfile' argument is missing. Check input.";
        $self->help and exit(0);
    }
    elsif ($opt->{all} && ! -e $opt->{gff}) {
        say STDERR "\nERROR: The '--gff' file does not appear to exist. Check input.";
        $self->help and exit(0);
    }
    elsif (! $opt->{indir} && ! $opt->{all}) {
        say STDERR "\nERROR: The '--indir' option must be given if no gff file and '--all' option is given. Check input.";
        $self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _calculate_tir_stats($opt);
}

sub _calculate_tir_stats {
    my ($opt) = @_;

    my %tirstats_opts = (
	genome  => $opt->{genome},
	outfile => $opt->{outfile},
    );
    
    $tirstats_opts{all}       = $opt->{all} // 0;
    $tirstats_opts{dir}       = $opt->{indir} // 0;
    $tirstats_opts{gff}       = $opt->{gff} // 0;
    $tirstats_opts{subs_rate} = $opt->{subs_rate} // 1e-8;
    $tirstats_opts{threads}   = $opt->{threads} // 1;
    $tirstats_opts{clean}     = $opt->{clean} // 0;

    my $stats_obj = Tephra::TIR::TIRStats->new(%tirstats_opts);

    $stats_obj->calculate_tir_ages;
}

sub help {
    print STDERR<<END

  USAGE: tephra tirage [-h] [-m]
      -m --man      :   Get the manual entry for a command.
      -h --help     :   Print the command usage.
    
  Required:
      -g|genome     :   The genome sequences in FASTA format used to search for TIRs.
      -f|gff        :   The GFF3 file of TIRs in <--genome>.
      -o|outfile    :   The output file containing the age of each element.
      -i|indir      :   The input directory of superfamily exemplars.

  Options:
      -r|subs_rate  :   The nucleotide substitution rate to use (Default: 1e-8).
      -t|threads    :   The number of threads to use for clustering coding domains (Default: 1).
      -c|clean      :   Clean up all the intermediate files from PAML and clustalw (Default: yes).
      -a|all        :   Calculate age of all TIRs in <gff> instead of exemplars in <indir>.   
   
END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra tirage - Calculate the age distribution of TIR transposons.

=head1 SYNOPSIS    

 tephra tirage -g ref.fas -f ref_tephra_mutator.gff3 -o ref_classified_tirs -t 12 --clean

=head1 DESCRIPTION

 This subcommand takes a GFF3 of TIR transposons as input from Tephra and calculates
 the insertion time and age of each element using a substitution rate and model of evolution.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for TIRs.

=item -f, --gff

 The GFF3 file of TIRs in <--genome> as output by the 'tephra classifytirs' command.

=item -o, --outdir

 The output directory for placing categorized elements.

=back

=head1 OPTIONS

=over 2

=item -r, --subs_rate

 The nucleotide substitution rate to use for calculating insertion times (Default: 1e-8).

=item -t, --threads

 The number of threads to use for clustering coding domains (Default: 1).

=item -c, --clean

 Clean up all the intermediate files from PAML and clustalw (Default: yes).

=item -a, --all

 Calculate age of all TIRs in <gff> instead of exemplars in <indir>.

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
