package Tephra::Command::classifyltrs;
# ABSTRACT: Classify LTR retrotransposons into superfamilies.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::Classify::LTRSfams;
use Tephra::Classify::LTRFams;
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use File::Path          qw(make_path remove_tree);
use File::Basename;
use File::Spec;

sub opt_spec {
    return (    
	[ "genome|g=s",   "The genome sequences in FASTA format to search for LTR-RTs "   ],
	[ "repeatdb|d=s", "The file of repeat sequences in FASTA format to use for classification " ], 
	[ "gff|f=s",      "The GFF3 file of LTR-RTs in <genome> "    ],
	[ "outdir|o=s",   "The output directory for placing categorized elements " ],
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
    elsif (!$opt->{genome} || !$opt->{repeatdb} || !$opt->{gff}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $some = _classify_ltr_predictions($opt);
}

sub _classify_ltr_predictions {
    my ($opt) = @_;

    my $genome   = $opt->{genome};
    my $repeatdb = $opt->{repeatdb};
    my $gff      = $opt->{gff};
    my $outdir   = $opt->{outdir};

    unless ( -d $outdir ) {
	make_path( $outdir, {verbose => 0, mode => 0771,} );
    }
    
    my $classify_obj = Tephra::Classify::LTRSfams->new( 
	genome   => $genome, 
	repeatdb => $repeatdb, 
	gff      => $gff 
    );

    my ($header, $features) = $classify_obj->collect_gff_features($gff);

    my ($gypsy, $copia) = $classify_obj->find_gypsy_copia($features);
    my ($unc_fas, $ltr_rregion_map) = $classify_obj->find_unclassified($features);

    my $blast_out = $classify_obj->search_unclassified($unc_fas);
    $classify_obj->annotate_unclassified($blast_out, $gypsy, $copia, $features, $ltr_rregion_map);
    my $gyp_gff = $classify_obj->write_gypsy($gypsy, $header);
    my $cop_gff = $classify_obj->write_copia($copia, $header);
    my $unc_gff = $classify_obj->write_unclassified($features, $header);

    my $classify_fams_obj = Tephra::Classify::LTRFams->new(
	genome   => $genome,
	outdir   => $outdir,
    );

    $classify_fams_obj->extract_features($gyp_gff);
    $classify_fams_obj->extract_features($cop_gff);
    $classify_fams_obj->extract_features($unc_gff);

    #my $vmatch_args = $classify_fams_obj->collect_feature_args;
    #dd $vmatch_args and exit;
    my $vmatch_clusters = $classify_fams_obj->cluster_features;
    #$classify_fams_obj->parse_clusters($vmatch_clusters);
}

sub help {
    print STDERR<<END

USAGE: tephra classifyltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for LTR-RTs. 
    -d|repeatdb   :   The file of repeat sequences in FASTA format to use for classification. 
    -f|gff        :   The GFF3 file of LTR-RTs in <--genome>.

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra classifyltrs - 

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
