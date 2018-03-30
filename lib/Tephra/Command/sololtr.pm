package Tephra::Command::sololtr;
# ABSTRACT: Find solo-LTRs in a genome assembly.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use File::Spec;
use File::Basename;
use Tephra -command;
use Tephra::Genome::SoloLTRSearch;

sub opt_spec {
    return (
	[ "indir|i=s",        "The directory of LTR families in FASTA format. "                             ],
	[ "genome|g=s",       "The genome sequences in FASTA format to search for solo-LTRs "               ],
	[ "percentid|p=i",    "Percent identity threshold for matches. (Default 39)."                       ], 
	[ "percentcov|c=i",   "Percent coverage threshold for matches. (Default 80)."                       ],
	[ "matchlen|l=i",     "Length threshold for matches. (Default 80)."                                 ],
	[ "outfile|o=s",      "The GFF file to save results to. "                                           ],
	[ "report|r=s",       "Parse hmmsearch of each sequence and produce a summary of align statistics." ],
	[ "numfamilies|n=i",  "The number of families to analyze (Default: the top 20)."                    ],
	[ "allfamilies|a",    "Analyze all LTR-RT families for solo-LTRs (Default: no)."                    ],
	[ "seq|s=s",          "Extract query sequence from domain alignment."                               ],
	[ "threads|t=i",      "The number of threads to use for clustering coding domains (Default: 1)"     ],
	[ "clean",            "Clean up the intermediate alignment files (Default: yes) "                   ],
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
	$self->help and exit(0);
    }
    elsif (!$opt->{indir} || !$opt->{genome} || !$opt->{outfile}) {
	say STDERR "\n[ERROR]: Required arguments not given.\n";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{indir}) {
	say STDERR "\n[ERROR]: The '--indir' directory does not appear to exist. Check input.\n";
        $self->help and exit(0);
    }
    elsif (! -e $opt->{genome}) {
	say STDERR "\n[ERROR]: The '--genome' file does not appear to exist. Check input.\n";
        $self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _calculate_soloLTR_abund($opt);
}

sub _calculate_soloLTR_abund {
    my ($opt) = @_;

    my ($name, $path, $suffix) = fileparse($opt->{genome}, qr/\.[^.]*/);
    my $report = $opt->{report} // File::Spec->catfile($path, $name.'_tephra_soloLTRs.tsv');

    my $pid     = $opt->{percentid} // 39;
    my $pcov    = $opt->{percentcov} // 80;
    my $len     = $opt->{matchlen} // 80;
    my $seq     = $opt->{seq} // 0;
    my $famnum  = $opt->{numfamilies} // 20;
    my $threads = $opt->{threads} // 1;
    my $clean   = $opt->{clean} // 1;
    my $all     = $opt->{allfamilies} // 0;

    my $ill_obj = Tephra::Genome::SoloLTRSearch->new(
	dir          => $opt->{indir},
	genome       => $opt->{genome},
	outfile      => $opt->{outfile},
	report       => $report,
	percentid    => $pid,
	matchlen     => $len,
	numfamilies  => $famnum,
	seqfile      => $seq,
	threads      => $threads,
	clean        => $clean,
	fullanalysis => $all,
    );
    
    $ill_obj->find_soloLTRs;
}

sub help {
    my $desc = capture_merged {
	pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
	    -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
 USAGE: tephra sololtr [-h] [-m]
     -m --man         :    Get the manual entry for a command.
     -h --help        :    Print the command usage.

 Required:
     -i|indir         :    The input directory of LTR families in FASTA format.
     -g|genome        :    Input FASTA file of sequences to search.
     -o|outfile       :    The GFF file to write the solo-LTRs to.

 Options:
     -p|percentid     :    Percent identity threshold for matches. (Default 39).
                           NB: For a threshold of 80% say 80.
     -c|percentcov    :    Percent coverage threshold for matches. (Default 80).
                           NB: DEPRECATED Option will be ignored in v0.07.0+.
     -l|matchlen      :    Length threshold for matches. (Default 80).
                           NB: For a threshold of 80 bp say 80.
     -r|report        :    Parse hmmsearch of each sequence and produce a summary of align statistics.
     -s|seq           :    Extract query sequence from domain alignment.
     -n|numfamilies   :    The number of families to analyze (Default: the top 20).
     -a|allfamilies   :    Analyze all LTR-RT families for solo-LTRs (Default: no).      
     -t|threads       :    The number of threads to use for clustering coding domains (Default: 1).
     --clean          :    Clean up the intermediate alignment and HMMER files (Default: yes).

END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra sololtr - Find solo LTRs in a genome assembly

=head1 SYNOPSIS    

 tephra sololtr -i ref_tephra_gypsy_dir -g ref_masked.fas -r sololtr_report.tsv -l 80 -p 39 -f 90 -s

=head1 DESCRIPTION

 The recombination of LTR retrotransposons leaves hallmark solo-LTR fragments in a genome. This
 command identifies the abundance and nature of solo LTRs in a genome arising from each LTR family.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

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

=item -p, --percentid

 Percent identity threshold for matches. (default 39).
 NB: For a threshold of 80 percent say 80.

=item -f, --percentcov

 Percent coverage threshold for matches. (default 0.80).
 NB: For a threshold of 80 percent say 80. This option has
 been deprecated.

=item -l, --matchlen

 Length threshold for matches. (default 80).
 NB: For a threshold of 80 bp say 80.

=item -r, --report

 Parse hmmsearch of each sequence and produce a summary of align statistics.

=item -s, --seq

 Extract query sequence from domain alignment.

=item -n, --numfamilies

 The number of families to analyze (Default: the top 20). For smaller genomes like Arabidopsis, it
 is fine to use the --allfamilies option because the analysis will complete in less than a day. For large
 genomes like sunflower this would take months, so the best approach, for now, is to restrict the
 analysis to the 20 largest families, which is the default. This option allows you to expand or reduce the 
 number of families to analyze.
 
=item -a, --allfamilies

 Analyze all LTR-RT families for solo-LTRs (Default: no). As mentioned above, this is fine for small genomes,
 but can lead to very long run times for plant genomes. My advice is to run the analysis under default conditions
 first, then expand the scope of the analysis if it seems feasible.

=item -t, --threads

 The number of threads to use for clustering coding domains (Default: 1).

=item --clean

 Clean up the intermediate alignment and HMMER files (Default: yes).

=item -h, --help

 Print a usage statement.

=item -m, --man

 Print the full documentation.

=back

=cut
