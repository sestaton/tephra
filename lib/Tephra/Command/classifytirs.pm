package Tephra::Command::classifytirs;
# ABSTRACT: Classify TIR transposons into superfamilies and families.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Cwd           qw(abs_path);
use File::Path    qw(make_path);
use File::Basename;
use Tephra -command;
use Tephra::Classify::TIRSfams;
use Tephra::Classify::Fams;
#use Log::Any qw($log);
#use Data::Dump::Color;

sub opt_spec {
    return (    
	[ "genome|g=s",     "The genome sequences in FASTA format used to search for TIRs "                          ],
        [ "logfile=s",      "The file to use for logging results in addition to the screen "                         ],
        [ "repeatdb|d=s",   "The file of repeat sequences in FASTA format to use for classification "                ], 
        [ "hitlength|l=i",  "The alignment length cutoff for BLAST hits to the repeat database (Default: 80) "       ],
        [ "percentid|p=i",  "The percent identity cutoff for BLAST hits to the repeat database (Default: 80) "       ],
        [ "outgff|o=s",     "The output GFF3 file of classified TIRs in <genome>    "                                ],
        [ "ingff|i=s",      "The input GFF3 file of TIRs in <genome> "                                               ],
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
        say STDERR "\nERROR: Required arguments not given.\n";
        $self->help and exit(0);
    }
    elsif (! -e $opt->{genome}) {
        say STDERR "\nERROR: The genome file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
    elsif (! -e $opt->{repeatdb}) {
        say STDERR "\nERROR: The repeat database file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
    elsif (! -e $opt->{ingff}) {
        say STDERR "\nERROR: The input GFF3 file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my ($gffs, $log) = _classify_tir_superfamilies($opt);
    _classify_tir_families($opt, $gffs, $log);
}

sub _classify_tir_superfamilies {
    my ($opt) = @_;
    
    my ($name, $path, $suffix) = fileparse($opt->{ingff}, qr/\.[^.]*/);
    my $sf_gff = File::Spec->catfile( abs_path($path), $name.'_superfamilies'.$suffix );

    my %classify_opts = (
	genome   => $opt->{genome},
        gff      => $opt->{ingff},
        outfile  => $sf_gff
    );

    $classify_opts{logfile} = $opt->{logfile} if $opt->{logfile};
    my $classify_obj = Tephra::Classify::TIRSfams->new(%classify_opts);
    
    my ($logfile, $log);
    if ($opt->{logfile}) {
        $log = $classify_obj->get_tephra_logger($opt->{logfile});
    }
    else {
        my ($name, $path, $suffix) = fileparse($opt->{genome}, qr/\.[^.]*/);
        $logfile = File::Spec->catfile( abs_path($path), $name.'_tephra_classifytirs.log' );
        $log = $classify_obj->get_tephra_logger($logfile);
        say STDERR "\nWARNING: '--logfile' option not given so results will be appended to: $logfile.";
    }

    my $index = $classify_obj->index_ref($opt->{genome});
    my ($header, $features) = $classify_obj->collect_gff_features($opt->{ingff});

    my $all_ct = (keys %$features);
    my ($tcmoutfile, $tcmfas, $tc1_ct) = $classify_obj->find_tc1_mariner($features, $header, $index, $log);
    my ($hatoutfile, $hatfas, $hat_ct) = $classify_obj->find_hat($features, $header, $index, $log);
    my ($mutoutfile, $mutfas, $mut_ct) = $classify_obj->find_mutator($features, $header, $index, $log);
    my ($cacoutfile, $cacfas, $cac_ct) = $classify_obj->find_cacta($features, $header, $index, $log);
    my ($uncoutfile, $uncfas, $unc_ct) = $classify_obj->write_unclassified_tirs($features, $header, $index, $log);

    my @fastas = grep { defined && /\.fasta$/ } ($tcmfas, $hatfas, $mutfas, $cacfas, $uncfas);
    my @gffs   = grep { defined && /\.gff3$/  } ($tcmoutfile, $hatoutfile, $mutoutfile, $cacoutfile, $uncoutfile);

    my %gffs;
    for my $file (@gffs) {
	if ($file =~ /tc1/i) {
	    $gffs{'mariner'} = $file;
	}
	if ($file =~ /hat/i) {
	     $gffs{'hat'} = $file;
	}
	if ($file =~ /mut/i) {
	     $gffs{'mutator'} = $file;
	}
	if ($file =~ /cac/i) {
	     $gffs{'cacta'} = $file;
	}
	if ($file =~ /unc/i) {
	     $gffs{'unclassified'} = $file;
	}
    }

    if (@fastas && @gffs) {
	#my %outfiles = (
	    #fastas => \@fastas,
	    #gffs   => \@gffs
	 #);

	#$classify_obj->write_combined_output(\%outfiles);
	
	$log->info("Results - Total number of TIR elements:                   $all_ct");
	$log->info("Results - Number of Tc1-Mariner elements:                 $tc1_ct");
	$log->info("Results - Number of hAT elements:                         $hat_ct");
	$log->info("Results - Number of Mutator elements:                     $mut_ct");
	$log->info("Results - Number of CACTA elements:                       $cac_ct");
	$log->info("Results - Number of unclassified TIR elements:            $unc_ct");
	unlink $_ for @fastas;

	return (\%gffs, $log);
    }
    else {
	say STDERR "\nWARNING: No TIR elements were classified. Check input.\n";
    }
}

sub _classify_tir_families {
    my ($opt, $gffs, $log) = @_;

    my $threads = $opt->{threads} // 1;
    my $hpcov   = $opt->{percentcov} // 50;
    my $hpid    = $opt->{percentid} // 80;
    my $hlen    = $opt->{hitlen} // 80;
    my $debug   = $opt->{debug} // 0;
    my $type    = 'TIR';

    unless ( -d $outdir ) {
        make_path( $outdir, {verbose => 0, mode => 0771,} );
    }

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
    $classify_fams_obj->annotate_gff($annot_ids, $ingff);
    
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
USAGE: tephra classifytirs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format used to search for TIRs. 
    -d|repeatdb   :   The file of repeat sequences in FASTA format to use for classification. 
    -o|outgff     :   The output GFF3 file of classified TIRs in <genome>.
    -i|ingff      :   The input GFF3 file of TIRs in <genome>.
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
                                                                       
 tephra classifytirs - Classify TIR transposons into superfamilies.

=head1 SYNOPSIS    

   tephra classifytirs -g ref.fas -d repeatdb.fas -i ref_tephra.gff3 -o ref_tephra_classified.gff3 -r ref_classified_tirs -t 12

=head1 DESCRIPTION
 
 This subcommand takes a GFF3 as input from Tephra and classifies the TIR TEs into
 superfamilies and families.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

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
