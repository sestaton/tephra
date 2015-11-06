package Tephra::Command::findtrims;
# ABSTRACT: Find TRIM retrotransposons in a genome assembly.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::TRIM::TRIMSearch;
use Tephra::LTR::LTRRefine;
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use File::Basename;
use File::Spec;

sub opt_spec {
    return (    
	[ "genome|g=s",  "The genome sequences in FASTA format to search for LTR-RTs "   ],
	[ "trnadb|t=s",  "The file of tRNA sequences in FASTA format to search for PBS " ], 
	[ "hmmdb|d=s",   "The HMM db in HMMERv3 format to search for coding domains "    ],
	[ "clean",       "Clean up the index files (Default: yes) "                      ],
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

    my ($relaxed_gff, $strict_gff) = _run_trim_search($opt);
    if ($relaxed_gff && $strict_gff) {
	my $some = _refine_trim_predictions($relaxed_gff, $strict_gff, $opt->{genome});
    }
}

sub _refine_trim_predictions {
    my ($relaxed_gff, $strict_gff, $fasta) = @_;

    my $refine_obj = Tephra::LTR::LTRRefine->new( genome  => $fasta );
	
    my $relaxed_features
	= $refine_obj->collect_features({ gff => $relaxed_gff, pid_threshold => 85 });
    my $strict_features
	= $refine_obj->collect_features({ gff => $strict_gff,  pid_threshold => 99 });

    my $best_elements = $refine_obj->get_overlaps({ relaxed_features => $relaxed_features, 
						    strict_features  => $strict_features });
    
    my $combined_features = $refine_obj->reduce_features({ relaxed_features => $relaxed_features, 
							   strict_features  => $strict_features,
							   best_elements    => $best_elements });

    $refine_obj->sort_features({ gff               => $relaxed_gff, 
				 combined_features => $combined_features });
}

sub _run_trim_search {
    my ($opt) = @_;
    
    my $genome = $opt->{genome};
    my $hmmdb  = $opt->{hmmdb};
    my $trnadb = $opt->{trnadb};
    my $clean  = defined $opt->{clean} ? $opt->{clean} : 0;
    
    my $trim_search = Tephra::TRIM::TRIMSearch->new( 
	genome => $genome, 
	hmmdb  => $hmmdb,
	trnadb => $trnadb, 
	clean  => $clean
    );

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $index = $genome.".index";

    my @suff_args = qq(-db $genome -indexname $index -tis -suf -lcp -ssp -sds -des -dna);
    $trim_search->create_index(\@suff_args);
    
    my $strict_gff  = $trim_search->trim_search_strict($index);
    my $relaxed_gff = $trim_search->trim_search_relaxed($index);

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
    -c|clean      :   Clean up the index files (Default: yes).

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

=item -t, --trnadb

 The file of tRNA sequences in FASTA format to search for PBS.

=item -d, --hmmdb

 The HMM db in HMMERv3 format to search for coding domains.

=back

=head1 OPTIONS

=over 2

=item -c, --clean

 Clean up the index files (Default: yes).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
