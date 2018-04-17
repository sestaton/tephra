package Tephra::Command::classifyltrs;
# ABSTRACT: Classify LTR retrotransposons into superfamilies and families.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Cwd           qw(abs_path);
use File::Path    qw(make_path remove_tree);
use File::Basename;
use Tephra -command;
use Tephra::Classify::LTRSfams;
use Tephra::Classify::Fams;

sub opt_spec {
    return (    
	[ "genome|g=s",     "The genome sequences in FASTA format used to search for LTR-RTs "                       ],
	[ "logfile=s",      "The file to use for logging results in addition to the screen "                         ],
	[ "repeatdb|d=s",   "The file of repeat sequences in FASTA format to use for classification "                ], 
	[ "hitlength|l=i",  "The alignment length cutoff for BLAST hits to the repeat database (Default: 80) "       ],
	[ "percentid|p=i",  "The percent identity cutoff for BLAST hits to the repeat database (Default: 80) "       ],
	[ "outgff|o=s",     "The output GFF3 file of classified LTR-RTs in <genome> "                                ],
	[ "ingff|i=s",      "The input GFF3 file of LTR-RTs in <genome> "                                            ],
	[ "outdir|r=s",     "The output directory for placing categorized elements "                                 ],
	[ "threads|t=i",    "The number of threads to use for clustering coding domains (Default: 1)  "              ],
	[ "percentcov|c=i", "The percent coverage cutoff for the shorter element in pairwise matches (Default: 50) " ],
        [ "percentid|p=i",  "The percent identity cutoff for classification of pairwise matches (Default: 80) "      ],
        [ "hitlen|l=i",     "The minimum length for classifying pairwise BLAST hits (Default: 80) "                  ],
	[ "debug",          "Show external command for debugging (Default: no) "                                     ],
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
    elsif (!$opt->{genome} || !$opt->{repeatdb} || !$opt->{ingff} || !$opt->{outgff} || !$opt->{outdir}) {
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
    elsif (! -e $opt->{ingff}) {
        say STDERR "\n[ERROR]: The input GFF3 file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my ($gffs, $lard_index, $log) = _classify_ltr_superfamilies($opt);
    my $somb = _classify_ltr_families($opt, $gffs, $lard_index, $log);
}

sub _classify_ltr_superfamilies {
    my ($opt) = @_;

    my $threads = $opt->{threads} // 1;

    unless ( -d $opt->{outdir} ) {
	make_path( $opt->{outdir}, {verbose => 0, mode => 0771,} );
    }
    
    my %classify_opts = (
	genome   => $opt->{genome},
        repeatdb => $opt->{repeatdb},
        gff      => $opt->{ingff},
        threads  => $threads,
    );

    $classify_opts{logfile} = $opt->{logfile} if $opt->{logfile};
    my $classify_obj = Tephra::Classify::LTRSfams->new(%classify_opts);

    my ($logfile, $log);
    if ($opt->{logfile}) {
        #$logfile = $self->logfile;
        $log = $classify_obj->get_tephra_logger($opt->{logfile});
    }
    else {
        my ($name, $path, $suffix) = fileparse($opt->{genome}, qr/\.[^.]*/);
        #my $lname = $self->is_trim ? 'tephra_findtrims.log' : 'tephra_findltrs.log';
        $logfile = File::Spec->catfile( abs_path($path), $name.'_tephra_classifyltrs.log' );
        $log = $classify_obj->get_tephra_logger($logfile);
        say STDERR "\n[WARNING]: '--logfile' option not given so results will be appended to: $logfile.";
    }

    my ($header, $features) = $classify_obj->collect_gff_features($opt->{ingff});
    my ($gypsy, $copia) = $classify_obj->find_gypsy_copia($features);

    my ($unc_fas, $ltr_rregion_map) = $classify_obj->find_unclassified($features);
    my $blast_out = $classify_obj->search_unclassified($unc_fas);
    $classify_obj->annotate_unclassified($blast_out, $gypsy, $copia, $features, $ltr_rregion_map);
    unlink $unc_fas;

    my ($gyp_gff, $cop_gff, $unc_gff, $lard_index, %gffs);
    if (%$gypsy) {
	$gyp_gff = $classify_obj->write_gypsy($gypsy, $header, $log);
	$gffs{'gypsy'} = $gyp_gff;
    }

    if (%$copia) {
        $cop_gff = $classify_obj->write_copia($copia, $header, $log);
	$gffs{'copia'} = $cop_gff;
    }

    if (%$features) {
        ($unc_gff, $lard_index) = $classify_obj->write_unclassified($features, $header, $log);
	$gffs{'unclassified'} = $unc_gff;
    }

    return (\%gffs, $lard_index, $log);
}

sub _classify_ltr_families {
    my ($opt, $gffs, $lard_index, $log) = @_;

    my $threads = $opt->{threads} // 1;
    my $hpcov   = $opt->{percentcov} // 50;
    my $hpid    = $opt->{percentid} // 80;
    my $hlen    = $opt->{hitlen} // 80;
    my $debug   = $opt->{debug} // 0;
    my $type    = 'LTR';

    my $classify_fams_obj = Tephra::Classify::Fams->new(
	genome        => $opt->{genome},
	gff           => $opt->{outgff},
	outdir        => $opt->{outdir},
	threads       => $threads,
	blast_hit_cov => $hpcov,
        blast_hit_pid => $hpid,
        blast_hit_len => $hlen,
	debug         => $debug,
	type          => $type,
    );

    my ($outfiles, $annot_ids) = $classify_fams_obj->make_families($gffs, $log);
    $classify_fams_obj->combine_families($outfiles);
    $classify_fams_obj->annotate_gff({ annotated_ids => $annot_ids, 
				       annotated_idx => $lard_index, 
				       input_gff     => $opt->{ingff}, 
				       te_type       => 'LTR' });
    
    unlink $_ for values %$gffs;
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra classifyltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format used to search for LTR-RTs. 
    -d|repeatdb   :   The file of repeat sequences in FASTA format to use for classification. 
    -o|outgff     :   The output GFF3 file of classified LTR-RTs in <genome>.
    -i|ingff      :   The input GFF3 file of LTR-RTs in <genome>.
    -r|outdir     :   The output directory for placing categorized elements.

Options:
    -t|threads    :   The number of threads to use for clustering coding domains (Default: 1).    
    -c|percentcov :   The percent coverage cutoff for the shorter element in pairwise matches (Default: 50).
    -p|percentid  :   The percent identity cutoff for classification of pairwise matches (Default: 80).
    -l|hitlen     :   The minimum length for classifying pairwise BLAST hits (Default: 80).
    --debug       :   Show external commands for debugging (Default: no).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra classifyltrs - Classify LTR retrotransposons into superfamilies and families

=head1 SYNOPSIS    

 tephra classifyltrs -g ref.fas -d repeatdb.fas -i ref_tephra.gff3 -o ref_tephra_classified.gff3 -r ref_classified_ltrs -t 12

=head1 DESCRIPTION

 This subcommand takes a GFF3 as input from Tephra and classifies the LTR-RTs first into
 superfamilies, then into families based on shared features.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -d, --repeatdb

 The file of repeat sequences in FASTA format to use for classification.

=item -o, --outgff

 The output GFF3 file of classified LTR-RTs in <genome>.

=item -i, --ingff

 The input GFF3 file of LTR-RTs in <genome>.

=item -r, --outdir

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

=item --debug

 Show external commands for debugging (Default: no).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
