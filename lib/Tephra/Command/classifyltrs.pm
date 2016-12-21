package Tephra::Command::classifyltrs;
# ABSTRACT: Classify LTR retrotransposons into superfamilies and families.

use 5.014;
use strict;
use warnings;
use File::Path qw(make_path remove_tree);
use Tephra -command;
use Tephra::Classify::LTRSfams;
use Tephra::Classify::LTRFams;
use Tephra::LTR::MakeExemplars;

sub opt_spec {
    return (    
	[ "genome|g=s",     "The genome sequences in FASTA format used to search for LTR-RTs "                       ],
	[ "repeatdb|d=s",   "The file of repeat sequences in FASTA format to use for classification "                ], 
	[ "hitlength|l=i",  "The alignment length cutoff for BLAST hits to the repeat database (Default: 80) "       ],
	[ "percentid|p=i",  "The percent identity cutoff for BLAST hits to the repeat database (Default: 80) "       ],
	[ "gff|f=s",        "The GFF3 file of LTR-RTs in <genome> "                                                  ],
	[ "outdir|o=s",     "The output directory for placing categorized elements "                                 ],
	[ "threads|t=i",    "The number of threads to use for clustering coding domains "                            ],
	[ "percentcov|c=i", "The percent coverage cutoff for the shorter element in pairwise matches (Default: 50) " ],
        [ "percentid|p=i",  "The percent identity cutoff for classification of pairwise matches (Default: 80) "      ],
        [ "hitlen|l=i",     "The minimum length for classifying pairwise BLAST hits (Default: 80) "                  ],
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

    my $gffs = _classify_ltr_superfamilies($opt);
    my $somb = _classify_ltr_families($opt, $gffs);
}

sub _classify_ltr_superfamilies {
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
    
    return ({ gypsy => $gyp_gff, copia => $cop_gff, unclassified => $unc_gff });
}

sub _classify_ltr_families {
    my ($opt, $gffs) = @_;

    my $genome   = $opt->{genome};
    my $repeatdb = $opt->{repeatdb};
    my $gff      = $opt->{gff};
    my $outdir   = $opt->{outdir};
    my $threads  = $opt->{threads} // 1;
    my $hpcov    = $opt->{percentcov} // 50;
    my $hpid     = $opt->{percentid} // 80;
    my $hlen     = $opt->{hitlen} // 80;

    my ($cop_gff, $gyp_gff, $unc_gff) = @{$gffs}{qw(copia gypsy unclassified)};

    for my $file ($cop_gff, $gyp_gff, $unc_gff) {
	unless (defined $file && -e $file) {
	    say "\nERROR: There was an error generating GFF3 for family level classification. Exiting.\n";
	}
    }

    my $classify_fams_obj = Tephra::Classify::LTRFams->new(
	genome        => $genome,
	outdir        => $outdir,
	threads       => $threads,
	blast_hit_cov => $hpcov,
        blast_hit_pid => $hpid,
        blast_hit_len => $hlen,
    );

    my ($outfiles, $annot_ids) = $classify_fams_obj->make_ltr_families($gffs);
    $classify_fams_obj->combine_families($outfiles);
    $classify_fams_obj->annotate_gff($annot_ids, $gff);
    unlink $gyp_gff, $cop_gff, $unc_gff;

    #say STDERR '=' x 50;
    #say STDERR join "\t", 'Gypsy_families', 'Gypsy_singletons', 'Copia_families', 'Copia_singletons', 
        #'Unclassified_families', 'Unclassified_singletons';
    #say STDERR join "\t", @{$gyp_ct}{qw(family_count singleton_count)}, @{$cop_ct}{qw(family_count singleton_count)},
        #@{$unc_ct}{qw(family_count singleton_count)};
    #say STDERR '=' x 50;
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
    -c|percentcov :   The percent coverage cutoff for the shorter element in pairwise matches (Default: 50).
    -p|percentid  :   The percent identity cutoff for classification of pairwise matches (Default: 80).
    -l|hitlen     :   The minimum length for classifying pairwise BLAST hits (Default: 80).

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

=item -c, --percentcov

 The percent coverage cutoff for the shorter element in pairwise matches (Default: 50). This option is used for family-level
 classifications. Elements will be added to existing families if this criteria is met, among others. This threshold has the 
 largest impact on family size, so use carefully. Increase for stringency and smaller families, decrease for reduced stringency 
 and larger families.

=item -l, --hitlength

 The alignment length cutoff for pairwise BLAST hits (Default: 80). This option is used for family-level                  
 classifications. It is recommended to not alter this value, as it has little impact. If necessary, it is best to alter this 
 while holding the percent coverage and percent identity unchanged.

=item -p, --percentid

 The percent identity cutoff for classification of pairwise matches (Default: 80). This option is used for family-level
 classifications. Increase for stringency and smaller families, decrease for reduced stringency               
 and larger families. It is recommended to leave this value unchanged and only change the percent coverage cutoff.

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
