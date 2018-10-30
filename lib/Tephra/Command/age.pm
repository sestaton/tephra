package Tephra::Command::age;
# ABSTRACT: Calculate the age distribution of LTR or TIR transposons.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Tephra -command;
use Tephra::Stats::Age;

sub opt_spec {
    return (
	[ "genome|g=s",    "The genome sequences in FASTA format used to search for LTRs/TIRs "         ],
	[ "gff|f=s",       "The GFF3 file of LTRs/TIRs in <genome> "                                    ],
	[ "outfile|o=s",   "The output file containing the age of each element "                        ],
	[ "subs_rate|r=f", "The nucleotide substitution rate to use (Default: 1e-8) "                   ],
	[ "threads|t=i",   "The number of threads to use for clustering coding domains "                ],
	[ "indir|i=s",     "The input directory of superfamily exemplars "                              ],
	[ "all|a",         "Calculate age of all LTRs/TIRs in <gff> instead of exemplars in <indir> "   ],
	[ "clean|c=i",     "Clean up all the intermediate files from PAML and clustalw (Default: yes) " ],
	[ "type=s",        "Type of transposon to calculate age for (must be 'ltr' or 'tir') "          ],
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
        $self->help and exit(0);
    }
    elsif (! $opt->{genome} || ! -e $opt->{genome}) {
        say STDERR "\n[ERROR]: The '--genome' argument was not given or the file does not exist. Check input.\n";
        $self->help and exit(0);
    }
    elsif (! $opt->{outfile}) {
	say STDERR "\n[ERROR]: The '--outfile' argument is missing. Check input.\n";
        $self->help and exit(0);
    }
    elsif ($opt->{all} && ! -e $opt->{gff}) {
        say STDERR "\n[ERROR]: The '--gff' file does not appear to exist. Check input.\n";
        $self->help and exit(0);
    }
    elsif (! $opt->{indir} && ! $opt->{all}) {
        say STDERR "\n[ERROR]: The '--indir' option must be given if no GFF3 file and '--all' option is given. Check input.\n";
        $self->help and exit(0);
    }
    elsif (! $opt->{type} || $opt->{type} !~ /ltr|tir/i) {
	say STDERR "\n[ERROR]: The '--type' argument is missing or does not match one of 'ltr' or 'tir' (case insensitive). Check input.\n";
        $self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _calculate_age_stats($opt);
}

sub _calculate_age_stats {
    my ($opt) = @_;

    my %agestats_opts = (
	genome  => $opt->{genome},
	outfile => $opt->{outfile},
    );
    
    $agestats_opts{all}       = $opt->{all} // 0;
    $agestats_opts{dir}       = $opt->{indir} // 0;
    $agestats_opts{gff}       = $opt->{gff} // 0;
    $agestats_opts{subs_rate} = $opt->{subs_rate} // 1e-8;
    $agestats_opts{threads}   = $opt->{threads} // 1;
    $agestats_opts{clean}     = $opt->{clean} // 1;
    $agestats_opts{type}      = lc($opt->{type});

    my $stats_obj = Tephra::Stats::Age->new(%agestats_opts);

    $stats_obj->calculate_ages;
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
  USAGE: tephra age [-h] [-m]
      -m --man      :   Get the manual entry for a command.
      -h --help     :   Print the command usage.
    
  Required:
      -g|genome     :   The genome sequences in FASTA format used to search for LTRs/TIRs.
      -f|gff        :   The GFF3 file of LTRs/TIRs in <--genome>.
      -o|outfile    :   The output file containing the age of each element.
      -i|indir      :   The input directory of superfamily exemplars.
      --type        :   Type of transposon to calculate age for (must be 'ltr' or 'tir').

  Options:
      -i|indir      :   The input directory of superfamily exemplars.
      -r|subs_rate  :   The nucleotide substitution rate to use (Default: 1e-8).
      -t|threads    :   The number of threads to use for clustering coding domains (Default: 1).
      -c|clean      :   Clean up all the intermediate files from PAML and clustalw (Default: yes).
      -a|all        :   Calculate age of all LTRs/TIRs in <gff> instead of exemplars in <indir>.   
   
END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra age - Calculate the age distribution of LTR or TIR transposons.

=head1 SYNOPSIS    

 tephra age -g ref.fas -f ref_tephra_tirs.gff3 -o ref_classified_tirs -t 12 --clean

=head1 DESCRIPTION

 This subcommand takes a GFF3 of LTR or TIR transposons as input from Tephra and calculates
 the insertion time and age of each element using a substitution rate and model of evolution.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTRs/TIRs.

=item -f, --gff

 The GFF3 file of LTRs/TIRs in <--genome> as output by the 'tephra classifytirs' or 'tephra classifyltrs' commands.

=item -o, --outdir

 The output directory for placing categorized elements.

=item --type

 Type of transposon to calculate age for (must be 'ltr' or 'tir').

=back

=head1 OPTIONS

=over 2

=item -o, --outdir

 The output directory for placing categorized elements.

=item -r, --subs_rate

 The nucleotide substitution rate to use for calculating insertion times (Default: 1e-8).

=item -t, --threads

 The number of threads to use for clustering coding domains (Default: 1).

=item -c, --clean

 Clean up all the intermediate files from PAML and clustalw (Default: yes).

=item -a, --all

 Calculate age of all LTRs/TIRs in <gff> instead of exemplars in <indir>.

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
