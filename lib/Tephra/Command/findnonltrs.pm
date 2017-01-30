package Tephra::Command::findnonltrs;
# ABSTRACT: Find non-LTR retrotransposons in a genome assembly.

use 5.014;
use strict;
use warnings;
use File::Spec;
use Tephra -command;
use Tephra::NonLTR::NonLTRSearch;
use Tephra::NonLTR::GFFWriter;

sub opt_spec {
    return (    
	[ "genome|g=s",  "The genome sequences in FASTA format to search for non-LTR-RTs " ],
	[ "pdir|d=s",    "The directory to search for HMMs (configured automatically) "    ],
	[ "outdir|d=s",  "The location to place the results "                              ],
	[ "gff|o=s",     "The GFF3 outfile to place the non-LTRs found in <genome> "       ],
	[ "verbose|v",   "Display progress for each chromosome (Default: no) "             ],
	[ "help|h",      "Display the usage menu and exit. "                               ],
        [ "man|m",       "Display the full manual. "                                       ],
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
    elsif (!$opt->{genome} || !$opt->{gff}) {
	say STDERR "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{genome}) {
        say STDERR "\nERROR: The genome file does not exist. Check arguments.";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my ($gff) = _run_nonltr_search($opt);
}

sub _run_nonltr_search {
    my ($opt) = @_;
    
    my $genome  = $opt->{genome};
    my $outdir  = $opt->{outdir};
    my $gff     = $opt->{gff};
    my $pdir    = $opt->{pdir} // $ENV{TEPHRA_DIR} // File::Spec->catdir($ENV{HOME}, '.tephra');
    my $verbose = $opt->{verbose} // 0;

    my $nonltr_obj = Tephra::NonLTR::NonLTRSearch->new(
	genome  => $genome,
	outdir  => $outdir,
	pdir    => $pdir,
	verbose => $verbose );
    
    my ($genomedir, $outputdir) = $nonltr_obj->find_nonltrs;

    my $gff_obj = Tephra::NonLTR::GFFWriter->new(
	genome   => $genome,
        fastadir => $genomedir,
	outdir   => $outputdir,
	gff      => $gff );

    $gff_obj->write_gff;
}

sub help {
    print STDERR<<END

USAGE: tephra findnonltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for non-LTR-RTs. 
    -o|gff        :   The GFF3 outfile to place the non-LTRs found in <genome>.

Options:
    -d|outdir     :   The location to place the results.
    -p|pdir       :   Location of the HMM models (Default: configured automatically).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findnonltrs - Find non-LTRs retrotransposons in a genome assembly.

=head1 SYNOPSIS    

 tephra findnonltrs -g genome.fas -d ref_nonltrs_results -o genome_nonltrs.gff3

=head1 DESCRIPTION
 
 Find non-LTR retrotransposons in a reference genome, classify them into known superfamilies, 
 and generate a GFF file showing their locations and properties.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

The genome sequences in FASTA format to search for non-LTR-RTs.

=item o, --gff

 The GFF3 outfile to place the non-LTRs found in <genome>.

=back

=head1 OPTIONS

=over 2

=item -o, --outdir

 The directory name to place the resulting GFF file (and run the analysis).
=item -d, --pdir

 The directory to search for HMMs. This should be configured automatically during installation and this option should only have to be used by developers.

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
