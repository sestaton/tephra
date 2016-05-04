package Tephra::Command::findtirs;
# ABSTRACT: Find TIR transposons in a genome assembly.

use 5.010;
use strict;
use warnings;
use File::Find;
use File::Basename;
use Tephra -command;
use Tephra::TIR::TIRSearch;

sub opt_spec {
    return (    
	[ "genome|g=s",  "The genome sequences in FASTA format to search for TIRs "   ],
	[ "hmmdb|d=s",   "The HMM db in HMMERv3 format to search for coding domains " ],
	[ "outfile|o=s", "The final combined and filtered GFF3 file of TIRs "         ],
	[ "index|i=s",   "The suffixerator index to use for the TIR search "          ],
	[ "clean",       "Clean up the index files (Default: yes) "                   ],
	[ "debug",       "Show external command for debugging (Default: no) "         ],
	[ "help|h",      "Display the usage menu and exit. "                          ],
        [ "man|m",       "Display the full manual. "                                  ],
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
    elsif (!$opt->{genome} || !$opt->{hmmdb}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $gff = _run_tir_search($opt);
}

sub _run_tir_search {
    my ($opt) = @_;
    
    my $genome  = $opt->{genome};
    my $hmmdb   = $opt->{hmmdb};
    my $index   = $opt->{index};
    my $outfile = $opt->{outfile};
    my $clean   = $opt->{clean} // 0;
    my $debug   = $opt->{debug} // 0;

    my @indexfiles;
    if (defined $index) {
        my ($name, $path, $suffix) = fileparse($index, qr/\.[^.]*/);
        my @files;
	for my $suf ('.des', '.lcp', '.llv', '.md5', '.prj', '.sds', '.suf')  {
            push @files, $index.$suf;
        }

	my $matchstr = join "|", @files;
        find( sub { push @indexfiles, $File::Find::name if -f and /$matchstr/ }, $path );
    }

    my $tir_search = Tephra::TIR::TIRSearch->new( 
	genome  => $genome, 
	outfile => $outfile,
	hmmdb   => $hmmdb,
	clean   => $clean,
	debug   => $debug,
    );

    unless (defined $index && @indexfiles == 7) {
	my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
	$index = $genome.'.index';

	my @suff_args = qq(-db $genome -indexname $index -tis -suf -lcp -des -ssp -dna -mirrored);
	$tir_search->create_index(\@suff_args);
    }
    
    my $gff = $tir_search->tir_search($index);

    return $gff;
}

sub help {
    print STDERR<<END

USAGE: tephra findtirs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for TIR TEs. 
    -d|hmmdb      :   The HMM db in HMMERv3 format to search for coding domains.

Options:
    -o|outfile    :   The final combined GFF3 file of TIRs.
    -i|index      :   The suffixerator index to use for the TIR search.
    --clean       :   Clean up the index files (Default: yes).
    --debug       :   Show external commands for debugging (Default: no).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findtirs - Find TIR transposons in a genome assembly.

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

=item -o, --outfile

 The name of a GFF3 file to place the filtered/reformatted TIR elements.

=item -i, --index

 The suffixerator index to use for the LTR search.

=item --clean

 Clean up the index files (Default: yes).

=item --debug

 Show external commands for debugging (Default: no).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
