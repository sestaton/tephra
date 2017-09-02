package Tephra::Command::findtirs;
# ABSTRACT: Find TIR transposons in a genome assembly.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Cwd           qw(getcwd);
use File::Copy    qw(copy);
use File::Find;
use File::Basename;
use Tephra -command;
use Tephra::Config::Exe;
use Tephra::TIR::TIRSearch;

sub opt_spec {
    return (    
	[ "genome|g=s",  "The genome sequences in FASTA format to search for TIRs "       ],
	[ "hmmdb|d=s",   "The HMM db in HMMERv3 format to search for coding domains "     ],
	[ "outfile|o=s", "The final combined and filtered GFF3 file of TIRs "             ],
	[ "index|i=s",   "The suffixerator index to use for the TIR search "              ],
	[ "threads|t=i", "The number of threads to use for the TIR search (Default: 1) "  ],
	[ "logfile|l=s", "The file to use for logging results in addition to the screen " ],
	[ "clean",       "Clean up the index files (Default: yes) "                       ],
	[ "debug",       "Show external command for debugging (Default: no) "             ],
	[ "help|h",      "Display the usage menu and exit. "                              ],
        [ "man|m",       "Display the full manual. "                                      ],
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
    elsif (!$opt->{genome}) {
	say STDERR "\nERROR: Required arguments not given.\n";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $gff = _run_tir_search($opt);
}

sub _run_tir_search {
    my ($opt) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($tephra_hmmdb) = @{$config}{qw(hmmdb)};
    my $hmmdb = $opt->{hmmdb} // $tephra_hmmdb;
    my $dir   = getcwd();

    my $genome  = $opt->{genome};
    my $index   = $opt->{index};
    my $outfile = $opt->{outfile};
    my $clean   = $opt->{clean} // 0;
    my $debug   = $opt->{debug} // 0;
    my $threads = $opt->{threads} // 1;
    my $logfile = $opt->{logfile} // File::Spec->catfile($dir, 'tephra_findtirs.log');
    #say $logfile;

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
	threads => $threads,
	logfile => $logfile,
    );

    unless (defined $index && @indexfiles == 7) {
	my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
	$index = $genome.'.index';

	my @suff_args = qq(-db $genome -indexname $index -tis -suf -lcp -des -ssp -dna -mirrored);
	$tir_search->create_index(\@suff_args, $logfile);
    }
    
    my $gff = $tir_search->tir_search($index);
    unlink $logfile unless -s $logfile;

    return $gff;
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra findtirs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for TIR TEs. 

Options:
    -o|outfile    :   The final combined GFF3 file of TIRs.
    -i|index      :   The suffixerator index to use for the TIR search.
    -d|hmmdb      :   The HMM db in HMMERv3 format to search for coding domains.    
    -t|threads    :   The number of threads to use for the TIR search (Default: 1).
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

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for TIR TEs.

=back

=head1 OPTIONS

=over 2

=item -o, --outfile

 The name of a GFF3 file to place the filtered/reformatted TIR elements.

=item -t, --threads

 The number of threads to use for the TIR search. Specifically,
 how many processors to use by the 'gt' program (Default: 1).

=item -i, --index

 The suffixerator index to use for the LTR search.

=item -d, --hmmdb

 The HMM db in HMMERv3 format to search for coding domains.

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
