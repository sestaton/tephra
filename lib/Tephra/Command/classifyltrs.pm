package Tephra::Command::classifyltrs;
# ABSTRACT: Classify LTR retrotransposons into superfamilies and families.

use 5.010;
use strict;
use warnings;
use File::Path qw(make_path remove_tree);
use Tephra -command;
use Tephra::Classify::LTRSfams;
use Tephra::Classify::LTRFams;
use Tephra::LTR::MakeExemplars;

sub opt_spec {
    return (    
	[ "genome|g=s",     "The genome sequences in FASTA format used to search for LTR-RTs "                 ],
	[ "repeatdb|d=s",   "The file of repeat sequences in FASTA format to use for classification "          ], 
	[ "hitlength|l=i",  "The alignment length cutoff for BLAST hits to the repeat database (Default: 80) " ],
	[ "percentid|p=i",  "The percent identity cutoff for BLAST hits to the repeat database (Default: 80) " ],
	[ "gff|f=s",        "The GFF3 file of LTR-RTs in <genome> "                                            ],
	[ "outdir|o=s",     "The output directory for placing categorized elements "                           ],
	[ "threads|t=i",    "The number of threads to use for clustering coding domains "                      ],
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
    elsif (!$opt->{genome} || !$opt->{repeatdb} || !$opt->{gff}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _classify_ltr_predictions($opt);
}

sub _classify_ltr_predictions {
    my ($opt) = @_;

    my $genome   = $opt->{genome};
    my $repeatdb = $opt->{repeatdb};
    my $gff      = $opt->{gff};
    my $outdir   = $opt->{outdir};
    my $threads  = $opt->{threads} // 1;

    unless ( -d $outdir ) {
	make_path( $outdir, {verbose => 0, mode => 0771,} );
    }
    
    my $classify_obj = Tephra::Classify::LTRSfams->new( 
	genome   => $genome, 
	repeatdb => $repeatdb, 
	gff      => $gff,
	threads  => $threads,
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
	threads  => $threads,
    );

    my $gyp_dir = $classify_fams_obj->extract_features($gyp_gff);
    my $gyp_clusters = $classify_fams_obj->cluster_features($gyp_dir);
    my ($gyp_fams, $gyp_ids, $gyp_ct) = $classify_fams_obj->parse_clusters($gyp_clusters);
    
    my $cop_dir = $classify_fams_obj->extract_features($cop_gff);
    my $cop_clusters = $classify_fams_obj->cluster_features($cop_dir);
    my ($cop_fams, $cop_ids, $cop_ct) = $classify_fams_obj->parse_clusters($cop_clusters);
    
    my $unc_dir = $classify_fams_obj->extract_features($unc_gff);
    my $unc_clusters = $classify_fams_obj->cluster_features($unc_dir);
    my ($unc_fams, $unc_ids, $unc_ct) = $classify_fams_obj->parse_clusters($unc_clusters);

    my (%outfiles, %annot_ids);
    @outfiles{keys %$_}  = values %$_ for ($gyp_fams, $cop_fams, $unc_fams);
    @annot_ids{keys %$_} = values %$_ for ($gyp_ids, $cop_ids, $unc_ids);
    $classify_fams_obj->combine_families(\%outfiles);
    $classify_fams_obj->annotate_gff(\%annot_ids, $gff);

    my $gyp_exm_obj = Tephra::LTR::MakeExemplars->new(
	genome => $genome,
	dir    => $gyp_dir,
	gff    => $gyp_gff
    );

    $gyp_exm_obj->make_exemplars;

    my $cop_exm_obj = Tephra::LTR::MakeExemplars->new(
	genome => $genome,
	dir    => $cop_dir,
	gff    => $cop_gff
    );

    $cop_exm_obj->make_exemplars;
    unlink $gyp_gff, $cop_gff, $unc_gff;

    say STDERR '=' x 50;
    say STDERR join "\t", 'Gypsy_families', 'Gypsy_singletons', 'Copia_families', 'Copia_singletons', 
        'Unclassified_families', 'Unclassified_singletons';
    say STDERR join "\t", @{$gyp_ct}{qw(family_count singleton_count)}, @{$cop_ct}{qw(family_count singleton_count)},
        @{$unc_ct}{qw(family_count singleton_count)};
    say STDERR '=' x 50;
}

sub help {
    print STDERR<<END

USAGE: tephra classifyltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format used to search for LTR-RTs. 
    -d|repeatdb   :   The file of repeat sequences in FASTA format to use for classification. 
    -f|gff        :   The GFF3 file of LTR-RTs in <--genome>.
    -o|outdir     :   The output directory for placing categorized elements.

Options:
    -t|threads    :   The number of threads to use for clustering coding domains (Default: 1).    
    -l|hitlength  :   The alignment length cutoff for BLAST hits to the repeat database (Default: 80).
    -p|percentid  :   The percent identity cutoff for BLAST hits to the repeat database (Default: 80).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra classifyltrs - Classify LTR retrotransposons into superfamilies and families

=head1 SYNOPSIS    

 tephra classifyltrs -g ref.fas -d repeatdb.fas -f ref_tephra.gff3 -o ref_classified_ltrs -t 12

=head1 DESCRIPTION

 This subcommand takes a GFF3 as input from Tephra and classifies the LTR-RTs first into
 superfamilies, then into families based on shared features.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -d, --repeatdb

 The file of repeat sequences in FASTA format to use for classification.

=item -f, --gff

 The GFF3 file of LTR-RTs in <--genome> as output by the 'tephra findltrs' command.

=item -o, --outdir

 The output directory for placing categorized elements.

=back

=head1 OPTIONS

=over 2

=item -t, --threads

 The number of threads to use for clustering coding domains (Default: 1).

=item -l, --hitlength

 The alignment length cutoff for BLAST hits to the repeat database (Default: 80). For samples where you do not have a repeat database from a closely related species, it may help to lower this value to increase the number of classified elements.

=item -p, --percentid

 The percent identity cutoff for BLAST hits to the repeat database (Default: 80). For samples where you do not have a repeat database from a closely related species, it may help to lower this value to increase the number of classified elements.

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
