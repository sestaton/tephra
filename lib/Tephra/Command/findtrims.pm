package Tephra::Command::findtrims;
# ABSTRACT: Find TRIM retrotransposons in a genome assembly.

use 5.014;
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Spec;
use File::Basename;
use Tephra -command;
use Tephra::Config::Exe;
use Tephra::TRIM::TRIMSearch;
use Tephra::LTR::LTRRefine;

sub opt_spec {
    return (    
	[ "genome|g=s",  "The genome sequences in FASTA format to search for TRIMs "      ],
	[ "trnadb|t=s",  "The file of tRNA sequences in FASTA format to search for PBS "  ], 
	[ "hmmdb|d=s",   "The HMM db in HMMERv3 format to search for coding domains "     ],
	[ "outfile|o=s", "The final combined and filtered GFF3 file of TRIMs "            ],
	[ "logfile=s",   "The file to use for logging results in addition to the screen " ],
	[ "clean",       "Clean up the index files (Default: No) "                        ],
	[ "debug",       "Show external command for debugging (Default: No) "             ],
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
        $self->help;
        exit(0);
    }
    elsif (!$opt->{genome} || !$opt->{outfile}) {
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

    my ($relaxed_gff, $strict_gff) = _run_trim_search($opt);
    if ($relaxed_gff && $strict_gff) {
	my $some = _refine_trim_predictions($relaxed_gff, $strict_gff, $opt->{genome}, $opt->{outfile}, $opt->{logfile});
    }
    elsif ($relaxed_gff && !$strict_gff) {
	_write_unrefined_trims($opt, $relaxed_gff);
    }
    elsif (!$relaxed_gff && !$strict_gff) {
	say STDERR "\nWARNING: No TRIMs were found, so there will be no output.\n";
	exit(1);
    }
}

sub _refine_trim_predictions {
    my ($relaxed_gff, $strict_gff, $genome, $outfile, $logfile) = @_;

    my %refine_opts = ( genome => $genome, outfile => $outfile, is_trim => 1 );
    $refine_opts{logfile} = $logfile if $logfile;
    my $refine_obj = Tephra::LTR::LTRRefine->new(%refine_opts);
	
    my $relaxed_features = $refine_obj->collect_features({ gff => $relaxed_gff, pid_threshold => 85 });
    my $strict_features  = $refine_obj->collect_features({ gff => $strict_gff,  pid_threshold => 99 });

    my $best_elements = $refine_obj->get_overlaps({ relaxed_features => $relaxed_features, 
						    strict_features  => $strict_features });
    
    my $combined_features = $refine_obj->reduce_features({ relaxed_features => $relaxed_features, 
							   strict_features  => $strict_features,
							   best_elements    => $best_elements });

    $refine_obj->sort_features({ gff               => $relaxed_gff, 
				 combined_features => $combined_features });

    unlink $relaxed_gff, $strict_gff;
}

sub _run_trim_search {
    my ($opt) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($tephra_hmmdb, $tephra_trnadb) = @{$config}{qw(hmmdb trnadb)};

    my $genome = $opt->{genome};
    my $hmmdb  = $opt->{hmmdb} // $tephra_hmmdb;
    my $trnadb = $opt->{trnadb} // $tephra_trnadb;
    my $clean  = $opt->{clean} // 0;
    my $debug  = $opt->{debug} // 0;

    my $trim_search = Tephra::TRIM::TRIMSearch->new( 
	genome => $genome, 
	hmmdb  => $hmmdb,
	trnadb => $trnadb, 
	clean  => $clean,
	debug  => $debug,
    );

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $index = File::Spec->catfile( abs_path($path), $name.$suffix.'.index');

    $genome = abs_path($genome);
    my @suff_args = qq(-db $genome -indexname $index -tis -suf -lcp -ssp -sds -des -dna);
    $trim_search->create_index(\@suff_args);
    
    #my $strict_gff  = $trim_search->trim_search_strict($index);
    #my $relaxed_gff = $trim_search->trim_search_relaxed($index);
    my $strict_gff  = $trim_search->trim_search({ index => $index, mode => 'strict'  });
    my $relaxed_gff = $trim_search->trim_search({ index => $index, mode => 'relaxed' });

    return ($relaxed_gff, $strict_gff);
}

sub _write_unrefined_trims {
    my ($opt, $relaxed_gff) = @_;

    my %refine_opts = ( genome => $opt->{genome}, outfile => $opt->{outfile}, is_trim => 1 );
    $refine_opts{logfile} = $opt->{logfile} if $opt->{logfile};
    my $refine_obj = Tephra::LTR::LTRRefine->new(%refine_opts);

    $refine_obj->sort_features({ gff               => $relaxed_gff,
                                 combined_features => undef });
}

sub help {
    print STDERR<<END

USAGE: tephra findtrims [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for LTR-RTs. 
    -o|outfile    :   The final combined and filtered GFF3 file of TRIMs. 

Options:
    --logfile     :   The file to use for logging results in addition to the screen.
    -t|trnadb     :   The file of tRNA sequences in FASTA format to search for PBS. 
    -d|hmmdb      :   The HMM db in HMMERv3 format to search for coding domains.
    -c|clean      :   Clean up the index files (Default: No).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findtrims - Find TRIM retrotransposons in a genome assembly.

=head1 SYNOPSIS    

 tephra findtrims -g ref.fas -t trnadb.fas -p te_models.hmm

=head1 DESCRIPTION

 Terminal Repeats In Minature (TRIMs) are abundant in many eukaryotic genomes and this command
 allows you to identify the nature and properties of these elements in a genome. By comparing these
 to autonomous LTR retrotransposons it may be possible to understand their origin and abundance.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for LTR-RTs.

=item -o, --outfile

 The final combined and filtered GFF3 file of TRIMs.

=back

=head1 OPTIONS

=over 2

=item --logfile

 The file to use for logging results in addition to the screen.

=item -t, --trnadb

 The file of tRNA sequences in FASTA format to search for PBS.

=item -d, --hmmdb

 The HMM db in HMMERv3 format to search for coding domains.

=item -c, --clean

 Clean up the index files (Default: No).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
