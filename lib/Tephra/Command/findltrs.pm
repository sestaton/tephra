package Tephra::Command::findltrs;
# ABSTRACT: Find LTR retrotransposons in a genome assembly.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::LTR::LTRSearch;
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
	[ "hmmdb|p=s",   "The HMM db in HMMERv3 format to search for coding domains "    ],
	[ "outfile|o=s", "The final combined and filtered GFF3 file of LTR-RTs "         ],
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
    elsif (!$opt->{genome} || !$opt->{outfile}) { # || !$opt->{hmmdb}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my ($relaxed_gff, $strict_gff) = _run_ltr_search($opt);
    my $some = _refine_ltr_predictions($relaxed_gff, $strict_gff, $opt->{genome}, $opt->{outfile});
}

sub _refine_ltr_predictions {
    my ($relaxed_gff, $strict_gff, $fasta, $outfile) = @_;

    my $refine_obj = Tephra::LTR::LTRRefine->new(
	genome  => $fasta,
	outfile => $outfile,
    );

    my ($all_feats, $all_stats, $intervals)
	= $refine_obj->collect_features($relaxed_gff, '85');
    my ($part_feats, $part_stats, $part_int) 
	= $refine_obj->collect_features($strict_gff, '99');

    my $best_elements = $refine_obj->get_overlaps($all_feats, $part_feats, $intervals);
    
    my $combined_features = $refine_obj->reduce_features($all_feats, $part_feats, $best_elements, 
					                 $all_stats, $part_stats);

    # and return something
    $refine_obj->sort_features($relaxed_gff, $combined_features, $outfile);

}

sub _run_ltr_search {
    my ($opt) = @_;
    
    my $genome = $opt->{genome};
    my $hmmdb  = $opt->{hmmdb};
    my $trnadb = $opt->{trnadb};
    my $clean  = $opt->{clean};
    $clean //= 0;
    
    #say "testing clean: $clean" and exit;
    
    my $ltr_search = Tephra::LTR::LTRSearch->new( 
	genome => $genome, 
	hmmdb  => $hmmdb,
	trnadb => $trnadb, 
	clean  => $clean );

    my $strict_gff  = $ltr_search->ltr_search_strict;
    my $relaxed_gff = $ltr_search->ltr_search_relaxed;

    #my $exit_value = 1;

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
    -p|hmmdb      :   The HMM db in HMMERv3 format to search for coding domains.

Options:
    -c|clean      :   Clean up the index files (Default: yes).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findltrs - 

=head1 SYNOPSIS    

 tephra findltrs -i .. -n

=head1 DESCRIPTION


=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -p, --paired

A file of interleaved, paired reads in FASTA format.

=item -u, --unpaired

A file of unpaired reads in FASTA format.

=back

=head1 OPTIONS

=over 2

=item -t, --treads

The number of threads to use with VelvetOptimiser (Default: 1).

=item -s, --hashs

The starting hash length for Velvet (Default: 59).

=item -e, --hashe

The ending hash length for Velvet (Default: 89).

=item -h, --help

Print a usage statement. 

=item -m, --man

Print the full documentation.

=back

=cut
