package Tephra::Command::maskref;
# ABSTRACT: Mask a reference genome with transposons.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Tephra -command;
use Tephra::Genome::MaskRef;

sub opt_spec {
    return (    
	[ "genome|g=s",    "The genome sequences in FASTA format to search for LTR-RTs "                        ],
	[ "outfile|o=s",   "The output filename to use for the masked genome "                                  ],
	[ "repeatdb|d=s",  "The database of repeat sequences to use for masking "                               ],
	[ "percentid|p=i",  "The percent identity cutoff for classification of pairwise matches (Default: 80) " ],
	[ "hitlength|l=i", "The alignment length cutoff for hits to the repeat database (Default: 70) "         ],
	[ "threads|t=i",   "The number of threads to use for masking (Default: 1) "                             ],
	[ "splitsize|s=i", "The chunk size to process at a time (Default: 50kb) "                               ],
	[ "overlap|v=i",   "The overlap between the chunks of the chromosomes (Default: 100 bp) "               ],
	[ "clean|c=i",     "Clean up the index files (Default: yes) "                                           ],
	[ "help|h",        "Display the usage menu and exit. "                                                  ],
        [ "man|m",         "Display the full manual. "                                                          ],
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
    elsif (!$opt->{genome} || !$opt->{repeatdb}) {
	say STDERR "\n[ERROR]: Required arguments not given.\n";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{genome}) {
	say STDERR "\n[ERROR]: The genome file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
    elsif (! -e $opt->{repeatdb}) {
	say STDERR "\n[ERROR]: The repeat database file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $gff = _run_masking($opt);
}

sub _run_masking {
    my ($opt) = @_;
    
    my $genome    = $opt->{genome};
    my $repeatdb  = $opt->{repeatdb};
    my $outfile   = $opt->{outfile};
    my $percentid = $opt->{percentid} // 80;
    my $hitlength = $opt->{hitlength} // 70;
    my $clean     = $opt->{clean} // 1;
    my $threads   = $opt->{threads} // 1;
    my $splitsize = $opt->{splitsize} // 5e4;
    my $overlap   = $opt->{overlap} // 100;

    my $mask_obj = Tephra::Genome::MaskRef->new( 
	genome    => $genome, 
	repeatdb  => $repeatdb,
	hit_pid   => $percentid,
	hitlength => $hitlength,
	outfile   => $outfile,
	clean     => $clean,
	threads   => $threads,
	splitsize => $splitsize,
	overlap   => $overlap,
    );
    
    $mask_obj->mask_reference;
    #my ($masked_ref, $report, $genome_length) = $mask_obj->mask_reference;
    #$mask_obj->get_masking_results($report, $genome_length);

}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra maskref [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for TIR TEs. 
    -d|repeatdb   :   The database of repeat sequences to use for masking.

Options:
    -o|outfile    :   The output filename to use for the masked genome.
    -p|percentid  :   The percent identity cutoff for classification of pairwise matches (Default: 80).
    -l|hitlength  :   The alignment length cutoff for hits to the repeat database (Default: 70).
    -s|splitsize  :   The chunk size to process at a time (Default: 50kb).
    -v|overlap    :   The overlap between the chunks of the chromosomes (Default: 100 bp).
    -t|threads    :   The number of threads to use for masking (Default: 1).
    -c|clean      :   Clean up the index files (Default: yes).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra maskref - Mask a reference genome with a custom repeat library.

=head1 SYNOPSIS    

 tephra maskref -g ref.fas -d repeatdb.fas -t 12

=head1 DESCRIPTION

 Mask a reference genome with one type of transposons to reduce false positives, and
 search time, in subsequent searches for other transposon types.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for TIR TEs.

=item -d, --repeatdb

 The database of repeat sequences to use for masking.

=back

=head1 OPTIONS

=over 2

=item -o, --outfile

 The output filename to use for the masked genome. If not given, the output will be named "input_masked.fas."
 For example, if the input is "seqs.fas" the output filename would be "seqs_masked.fas."

=item -t, --threads

 The number of threads to use for masking (Default: 1). 

=item -p, --percentid

 The percent identity cutoff for classification of pairwise matches with the repeat database
 and the reference (Default: 80).

=item -l, --hitlength

 The alignment length cutoff for hits to the repeat database (Default: 70). Lower this value to increase the
 percent of the reference that is masked, or alternatively, increase to be more conservative with masking.

=item -s, --splitsize

 The chunk size to process at a time (Default: 50kb). Increasing the size to processes will speed up the
 execution time but will increase memory usage. 

=item -v, -- overlap 

 The overlap between the chunks of the chromosomes (Default: 100 bp). Increasing this value will slow down 
 processing. The goal is to reduce artifacts by spliting up the chromosomes.

=item -c, --clean

 Clean up the index files (Default: yes).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
