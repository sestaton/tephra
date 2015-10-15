package Tephra::Command::findnonltrs;
# ABSTRACT: Find non-LTR retrotransposons in a genome assembly.

use 5.010;
use strict;
use warnings;
use File::Basename;
use Tephra -command;
use Tephra::NonLTR::NonLTRSearch;
#use Tephra::LTR::LTRRefine;

sub opt_spec {
    return (    
	[ "genome|g=s",  "The genome sequences in FASTA format to search for non-LTR-RTs "    ],
	[ "hmmdir|d=s",  "The directory of HMMs in HMMERv3 format to search for coding domains "     ],
	[ "outdir|o=s",  "The location to place the results "          ],
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
    elsif (!$opt->{genome}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my ($relaxed_gff, $strict_gff) = _run_nonltr_search($opt);
    #my $some = _refine_ltr_predictions($relaxed_gff, $strict_gff, $opt);
}

sub _run_nonltr_search {
    my ($opt) = @_;
    
    my $genome = $opt->{genome};
    
    my $ltr_search = Tephra::NonLTR::NonLTRSearch->new( 
	genome => $genome, 
	hmmdb  => $hmmdb,
	trnadb => $trnadb, 
	clean  => $clean 
    );

    unless (defined $index) {
	my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
	$index = $genome.".index";
    
	my @suff_args = qq(-db $genome -indexname $index -tis -suf -lcp -ssp -sds -des -dna);
	$ltr_search->create_index(\@suff_args);
    }
    
    my $strict_gff  = $ltr_search->ltr_search_strict($index);
    my $relaxed_gff = $ltr_search->ltr_search_relaxed($index);

    return ($relaxed_gff, $strict_gff);
}

sub help {
    print STDERR<<END

USAGE: tephra findltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for LTR-RTs. 
    -t|trnadb     :   The file of tRNA sequences in FASTA format to search for PBS. 
    -d|hmmdb      :   The HMM db in HMMERv3 format to search for coding domains.

Options:
    -o|outfile    :   The final combined and filtered GFF3 file of LTR-RTs.
    -i|index      :   The suffixerator index to use for the LTR search.
    -r|dedup      :   Discard elements with duplicate coding domains (Default: no).
    --tnpfilter   :   Discard elements containing transposase domains (Default: no).
    -c|clean      :   Clean up the index files (Default: yes).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findltrs - Find LTR retrotransposons in a genome assembly.

=head1 SYNOPSIS    

 tephra findltrs -g ref.fas -t trnadb.fas -d te_models.hmm --tnpfilter --clean

=head1 DESCRIPTION
 
 Search a reference genome and find LTR-RTs under relaxed and strict conditions (more on
 this later...), filter all predictions by a number of criteria, and generate a consensus
 set of the best scoring elements.

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
