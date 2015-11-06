package Tephra::Command::findnonltrs;
# ABSTRACT: Find non-LTR retrotransposons in a genome assembly.

use 5.010;
use strict;
use warnings;
use File::Basename;
use Tephra -command;
use Tephra::NonLTR::NonLTRSearch;
use Tephra::NonLTR::GFFWriter;

sub opt_spec {
    return (    
	[ "genomedir|g=s", "The genome sequences in FASTA format to search for non-LTR-RTs " ],
	[ "pdir|d=s",      "The directory to search for HMMs (configured automatically) "    ],
	[ "outdir|o=s",    "The location to place the results "                              ],
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
    elsif (!$opt->{genomedir} || !$opt->{outdir}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my ($gff) = _run_nonltr_search($opt);
}

sub _run_nonltr_search {
    my ($opt) = @_;
    
    my $genomedir = $opt->{genomedir};
    my $outdir    = $opt->{outdir};
    my $pdir      = $opt->{pdir};
    $pdir //= File::Spec->catdir($ENV{HOME}, '.tephra');

    my $nonltr_obj = Tephra::NonLTR::NonLTRSearch->new(
	fastadir => $genomedir,
	outdir   => $outdir,
	pdir     => $pdir );

    $nonltr_obj->find_nonltrs;

    my $gff_obj = Tephra::NonLTR::GFFWriter->new(
        fastadir => $genomedir,
	outdir   => $outdir );

    $gff_obj->write_gff;
}

sub help {
    print STDERR<<END

USAGE: tephra findnonltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genomedir  :   The genome sequences in FASTA format (one sequence per file) 
                      to search for non-LTR-RTs. 
    -o|outdir     :   The location to place the results. 

Options:
    -p|pdir       :   Location of the HMM models (Default: configured automatically).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findnonltrs - Find non-LTRs retrotransposons in a genome assembly.

=head1 SYNOPSIS    

 tephra findnonltrs

=head1 DESCRIPTION
 


=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -t, --trnadb

 The file of tRNA sequences in FASTA format to search for PBS.

=item -d, --hmmdb

 The HMM db in HMMERv3 format to search for coding domains.

=back

=head1 OPTIONS

=over 2

=item -o, --outfile

 The final combined and filtered GFF3 file of LTR-RTs.

=item -i, --index

 The suffixerator index to use for the LTR search.

=item -r, --dedup

 Discard elements with duplicate coding domains (Default: no).

=item --tnpfilter

 Discard elements containing transposase domains (Default: no).

=item -c, --clean

 Clean up the index files (Default: yes).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
