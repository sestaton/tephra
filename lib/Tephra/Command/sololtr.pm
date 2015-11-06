package Tephra::Command::sololtr;
# ABSTRACT: Find solo-LTRs in a genome assembly.

use 5.010;
use strict;
use warnings;
use Tephra::Genome::SoloLTRSearch;
use Tephra -command;

sub opt_spec {
    return (
	[ "indir|i=s",        "The directory of LTR families in FASTA format. "                             ],
	[ "genome|g=s",       "The genome sequences in FASTA format to search for solo-LTRs "               ],
	[ "percentident|p=f", "Percent identity threshold for matches. (default 0.39)."                     ], 
	[ "percentcov|c=f",   "Percent coverage threshold for matches. (default 0.80)."                     ],
	[ "matchlen|l=f",     "Length threshold for matches. (default 80)."                                 ],
	[ "report|r=s",       "Parse hmmsearch of each sequence and produce a summary of align statistics." ],
	[ "seq|s",            "Extract query sequence from domain alignment."                               ],
	[ "clean",            "Clean up the intermediate alignment files (Default: yes) "                   ],
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
    elsif (!$opt->{indir} || !$opt->{genome}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $some = _calculate_soloLTR_abund($opt);
}

sub _calculate_soloLTR_abund {
    my ($opt) = @_;

    my $dir    = $opt->{indir};
    my $genome = $opt->{genome};
    my $report = $opt->{report};
    my $pid    = defined $opt->{percentident} ? $opt->{percentident} : 0.39;
    my $pcov   = defined $opt->{percentcov} ? $opt->{percentcov} : 0.80;
    my $len    = defined $opt->{matchlen} ? $opt->{matchlen} : 80;
    my $seq    = defined $opt->{seq} ? 1 : 0;
    my $clean  = defined $opt->{clean} ? 1 : 0;

    my $ill_obj = Tephra::Genome::SoloLTRSearch->new(
	dir          => $dir,
	genome       => $genome,
	report       => $report,
	percentident => $pid,
	percentcov   => $pcov,
	matchlen     => $len,
	seq          => $seq,
	clean        => $clean,
    );
    
    $ill_obj->find_soloLTRs;
}

sub help {
    print STDERR<<END

 USAGE: tephra sololtr [-h] [-m]
     -m --man         :    Get the manual entry for a command.
     -h --help        :    Print the command usage.

 Required:
     -i|indir         :    The input directory of LTR families in FASTA format.
     -g|genome        :    Input FASTA file of sequences to search.

 Options:
     -p|percentident  :    Percent identity threshold for matches. (default 0.39).
                           NB: For a threshold of 80 percent say 0.80.
     -r|report        :    Parse hmmsearch of each sequence and produce a summary of align statistics.
     -s|seq           :    Extract query sequence from domain alignment.

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

=back

=head1 OPTIONS

=over 2

=item -p, --percentident

 Percent identity threshold for matches. (default 0.39).
 NB: For a threshold of 80 percent say 0.80.
 
=item -r, --report

 Parse hmmsearch of each sequence and produce a summary of align statistics.

=item -s, --seq

 Extract query sequence from domain alignment.

=item -h, --help

 Print a usage statement.

=item -m, --man

 Print the full documentation.

=back

=cut
