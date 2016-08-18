package Tephra::Command::sololtr;
# ABSTRACT: Find solo-LTRs in a genome assembly.

use 5.010;
use strict;
use warnings;
#use Pod::Find     qw(pod_where);
#use Pod::Usage    qw(pod2usage);
#use Capture::Tiny qw(:all);
use Tephra -command;
use Tephra::Genome::SoloLTRSearch;

sub opt_spec {
    return (
	[ "indir|i=s",        "The directory of LTR families in FASTA format. "                             ],
	[ "genome|g=s",       "The genome sequences in FASTA format to search for solo-LTRs "               ],
	[ "percentident|p=f", "Percent identity threshold for matches. (default 0.39)."                     ], 
	[ "percentcov|c=f",   "Percent coverage threshold for matches. (default 0.80)."                     ],
	[ "matchlen|l=f",     "Length threshold for matches. (default 80)."                                 ],
	[ "outfile|o=s",      "The GFF file to save results to. "                                           ],
	[ "report|r=s",       "Parse hmmsearch of each sequence and produce a summary of align statistics." ],
	[ "seq|s",            "Extract query sequence from domain alignment."                               ],
	[ "clean",            "Clean up the intermediate alignment files (Default: yes) "                   ],
	[ "threads|t=i",      "The number of threads to use for clustering coding domains "                 ],
	[ "help|h",           "Display the usage menu and exit. "                                           ],
	[ "man|m",            "Display the full manual. "                                                   ],
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
    elsif (!$opt->{indir} || !$opt->{genome} || !$opt->{outfile}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{indir}) {
	say "\nERROR: The '--indir' directory does not appear to exist. Check input.";
        $self->help and exit(0);
    }
    elsif (! -e $opt->{genome}) {
	say "\nERROR: The '--genome' file does not appear to exist. Check input.";
        $self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _calculate_soloLTR_abund($opt);
}

sub _calculate_soloLTR_abund {
    my ($opt) = @_;

    my $dir     = $opt->{indir};
    my $genome  = $opt->{genome};
    my $report  = $opt->{report};
    my $outfile = $opt->{outfile};
    my $pid     = $opt->{percentident} // 0.39;
    my $pcov    = $opt->{percentcov} // 0.80;
    my $len     = $opt->{matchlen} // 80;
    my $seq     = $opt->{seq} // 0;
    my $threads = $opt->{threads} // 1;
    my $clean   = $opt->{clean} // 1;

    my $ill_obj = Tephra::Genome::SoloLTRSearch->new(
	dir          => $dir,
	genome       => $genome,
	report       => $report,
	outfile      => $outfile,
	percentident => $pid,
	percentcov   => $pcov,
	matchlen     => $len,
	seq          => $seq,
	threads      => $threads,
	clean        => $clean,
    );
    
    $ill_obj->find_soloLTRs;
}

sub help {
    #my $stdout = capture_merged {
	#pod2usage(-verbose => 99, -sections => "NAME|SYNOPSIS|DESCRIPTION", -exitval => "noexit",
	#    -input => pod_where({-inc => 1}, __PACKAGE__));
    #};
    #chomp $stdout;

    print STDERR<<END

 USAGE: tephra sololtr [-h] [-m]
     -m --man         :    Get the manual entry for a command.
     -h --help        :    Print the command usage.

 Required:
     -i|indir         :    The input directory of LTR families in FASTA format.
     -g|genome        :    Input FASTA file of sequences to search.
     -o|outfile       :    The GFF file to write the solo-LTRs to.

 Options:
     -p|percentident  :    Percent identity threshold for matches. (default 0.39).
                           NB: For a threshold of 80 percent say 0.80.
     -f|percentcov    :    Percent coverage threshold for matches. (default 0.80).
                           NB: For a threshold of 80 percent say 0.80.
     -l|matchlen      :    Length threshold for matches. (default 80).
                           NB: For a threshold of 80 percent say 0.80.
     -r|report        :    Parse hmmsearch of each sequence and produce a summary of align statistics.
     -s|seq           :    Extract query sequence from domain alignment.
     -t|threads       :    The number of threads to use for clustering coding domains.
     --clean          :    Clean up the intermediate alignment and HMMER files (Default: yes).

END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra sololtr - Find solo LTRs in a genome assembly

=head1 SYNOPSIS    

 tephra sololtr sololtr -i ref_tephra_gypsy_dir -g ref_masked.fas -r sololtr_report.tsv -l 80 -p -c 0.09 -s

=head1 DESCRIPTION

 The recombination of LTR retrotransposons leaves hallmark solo-LTR fragments in a genome. This
 command identifies the abundance and nature of solo LTRs in a genome arising from each LTR family.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -i, --indir

 The input directory of LTR families in FASTA format.

=item -g, --genome

 Input (masked) FASTA file of genomic sequences to search.

=item -o, --outfile

 The name of the GFF file to write the solo-LTR locations to.

=back

=head1 OPTIONS

=over 2

=item -p, --percentident

 Percent identity threshold for matches. (default 0.39).
 NB: For a threshold of 80 percent say 0.80.

=item -f, --percentcov

 Percent coverage threshold for matches. (default 0.80).
 NB: For a threshold of 80 percent say 0.80.

=item -l, --matchlen

 Length threshold for matches. (default 80).
 NB: For a threshold of 80 percent say 0.80.

=item -r, --report

 Parse hmmsearch of each sequence and produce a summary of align statistics.

=item -s, --seq

 Extract query sequence from domain alignment.

=item -t, --threads

 The number of threads to use for clustering coding domains.

=item --clean

 Clean up the intermediate alignment and HMMER files (Default: yes).

=item -h, --help

 Print a usage statement.

=item -m, --man

 Print the full documentation.

=back

=cut
