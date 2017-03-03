package Tephra::Command::ltrage;
# ABSTRACT: Calculate the age distribution of LTR retrotransposons.

use 5.014;
use strict;
use warnings;
use Tephra -command;
use Tephra::LTR::LTRStats;
#use Data::Dump::Color;

sub opt_spec {
    return (
	[ "genome|g=s",    "The genome sequences in FASTA format used to search for LTR-RTs "           ],
	[ "gff|f=s",       "The GFF3 file of LTR-RTs in <genome> "                                      ],
	[ "outfile|o=s",   "The output file containing the age of each element "                        ],
	[ "subs_rate|r=f", "The nucleotide substitution rate to use (Default: 1e-8) "                   ],
	[ "threads|t=i",   "The number of threads to use for clustering coding domains "                ],
	[ "indir|i=s",     "The input directory of classifed LTR elements "                             ],
	[ "all|a",         "Calculate age of all LTR-RTs in <gff> instead of exemplars in <indir> "     ],
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

    my $some = _calculate_ltr_stats($opt);
}

sub _calculate_ltr_stats {
    my ($opt) = @_;

    my %ltrstats_opts = (
	genome  => $opt->{genome},
	outfile => $opt->{outfile},
    );
    
    $ltrstats_opts{all}       = $opt->{all} // 0;
    $ltrstats_opts{dir}       = $opt->{indir} // 0;
    $ltrstats_opts{gff}       = $opt->{gff} // 0;
    $ltrstats_opts{subs_rate} = $opt->{subs_rate} // 1e-8;
    $ltrstats_opts{threads}   = $opt->{threads} // 1;
    $ltrstats_opts{clean}     = $opt->{clean} // 0;

    my $stats_obj = Tephra::LTR::LTRStats->new(%ltrstats_opts);

    $stats_obj->calculate_ltr_ages;
}

sub help {
    print STDERR<<END

  USAGE: tephra ltrage [-h] [-m]
      -m --man      :   Get the manual entry for a command.
      -h --help     :   Print the command usage.
    
  Required:
      -g|genome     :   The genome sequences in FASTA format used to search for LTR-RTs.
      -f|gff        :   The GFF3 file of LTR-RTs in <--genome>.
      -o|outfile    :   The output file containing the age of each element.

  Options:
      -i|indir      :   The input directory of classified LTR elements.
      -r|subs_rate  :   The nucleotide substitution rate to use (Default: 1e-8).
      -t|threads    :   The number of threads to use for clustering coding domains (Default: 1).
      -c|clean      :   Clean up all the intermediate files from PAML and clustalw (Default: yes).
      -a|all        :   Calculate age of all LTR-RTs in <gff> instead of exemplars in <indir>.   
   
END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra ltrage - Calculate the age distribution of LTR retrotransposons.

=head1 SYNOPSIS    

 tephra ltrage -g ref.fas -f ref_tephra_gypsy.gff3 -o ref_classified_ltrs -t 12 --clean --all

=head1 DESCRIPTION

 This subcommand takes a GFF3 of LTR retrotransposons as input from Tephra and calculates
 the insertion time and age of each element using a substitution rate and model of evolution.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -f, --gff

 The GFF3 file of LTR-RTs in <--genome> as output by the 'tephra classifyltrs' command.

=item -o, --outdir

 The output directory for placing categorized elements.

=back

=head1 OPTIONS

=over 2

=item -i|indir

 The input directory of classified LTR elements.

=item -r, --subs_rate

 The nucleotide substitution rate to use for calculating insertion times (Default: 1e-8).

=item -t, --threads

 The number of threads to use for clustering coding domains (Default: 1).

=item -c, --clean

 Clean up all the intermediate files from PAML and clustalw (Default: yes).

=item -a, --all

 Calculate age of all LTR-RTs in <gff> instead of exemplars in <indir>.

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
