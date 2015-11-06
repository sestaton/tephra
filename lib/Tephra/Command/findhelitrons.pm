package Tephra::Command::findhelitrons;
# ABSTRACT: Find Helitons in a genome assembly.

use 5.010;
use strict;
use warnings;
use File::Basename;
use Tephra -command;
use Tephra::Hel::HelSearch;

sub opt_spec {
    return (    
	[ "genome|g=s",              "The genome sequences in FASTA format to search for Helitrons "   ],
	[ "helitronscanner_dir|d=s", "The HelitronScanner directory "    ],
	[ "outfile|o=s",             "The final combined and filtered GFF3 file of Helitrons "         ],
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
    elsif (!$opt->{genome} || !$opt->{outfile}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $gff = _run_helitron_search($opt);
}

sub _run_helitron_search {
    my ($opt) = @_;
    
    my $genome     = $opt->{genome};
    my $hscan_dir  = $opt->{helitronscanner_dir};
    my $gff        = $opt->{outfile};
    $hscan_dir //= File::Spec->catfile($ENV{HOME}, '.tephra', 'helitronscanner', 'HelitronScanner');

    my $hel_search = Tephra::Hel::HelSearch->new( 
	genome  => $genome, 
	helitronscanner_dir => $hscan_dir,
	outfile => $gff
    );

    my $hel_seqs = $hel_search->find_helitrons;
    $hel_search->make_hscan_gff($hel_seqs);
    
    return $gff;
}

sub help {
    print STDERR<<END

USAGE: tephra findhelitrons [-h] [-m]
    -m --man                :   Get the manual entry for a command.
    -h --help               :   Print the command usage.

Required:
    -g|genome               :   The genome sequences in FASTA format to search for Helitrons.. 
    -o|outfile              :   The final combined and filtered GFF3 file of Helitrons.

Options:
    -d|helitronscanner_dir  :   The HelitronScanner directory containing the ".jar" files and Training Set.
                                This should be configured automatically upon a successful install.

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findhelitrons - Find Helitrons in a genome assembly.

=head1 SYNOPSIS    

 tephra findtirs -g ref.fas -d te_models.hmm

=head1 DESCRIPTION

 Find TIR transposons in a reference genome assembly.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for TIR TEs.

=item -d, --hmmdb

 The HMM db in HMMERv3 format to search for coding domains.

=back

=head1 OPTIONS

=over 2

=item -i, --index

 The suffixerator index to use for the LTR search.

=item -c, --clean

 Clean up the index files (Default: yes).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
