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
    elsif (-e $opt->{outdir}) {
        say STDERR "\n[ERROR]: The output directory exists. This may cause problems. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my ($gffs, $mite_index, $log) = _classify_tir_superfamilies($opt);
    _classify_tir_families($opt, $gffs, $mite_index, $log);
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
        say STDERR "\n[WARNING]: '--logfile' option not given so results will be appended to: $logfile.";
    }

    my $index = $classify_obj->index_ref($opt->{genome});
    my ($header, $features) = $classify_obj->collect_gff_features($opt->{ingff});

    my $all_ct = (keys %$features);
    my ($tcmoutfile, $tcmfas, $tc1_ct) = $classify_obj->find_tc1_mariner($features, $header, $index, $log);
    my ($hatoutfile, $hatfas, $hat_ct) = $classify_obj->find_hat($features, $header, $index, $log);
    my ($mutoutfile, $mutfas, $mut_ct) = $classify_obj->find_mutator($features, $header, $index, $log);
    my ($cacoutfile, $cacfas, $cac_ct) = $classify_obj->find_cacta($features, $header, $index, $log);
    my $unc_obj = $classify_obj->write_unclassified_tirs($features, $header, $index, $log);

    my @fastas = grep { defined && /\.fasta$/ } 
        ($tcmfas, $hatfas, $mutfas, $cacfas, $unc_obj->{unc_fasta}, $unc_obj->{mite_fasta});
    my @gffs = grep { defined && /\.gff3$/ } 
        ($tcmoutfile, $hatoutfile, $mutoutfile, $cacoutfile, $unc_obj->{unc_outfile}, $unc_obj->{mite_outfile});

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
	if ($file =~ /mite/i) {
	    $gffs{'mite'} = $file;
        }
    }

    if (@fastas && @gffs) {
	unlink $_ for @fastas;
	
	$unc_obj->{mite_count} //= 0;
	$unc_obj->{unc_count} //= 0;
	my $tot_str = sprintf("%-70s %-10s", "Results - Total number of TIR elements:", $all_ct);
	my $tc1_str = sprintf("%-70s %-10s", "Results - Number of Tc1-Mariner elements:", $tc1_ct);
	my $hat_str = sprintf("%-70s %-10s", "Results - Number of hAT elements:", $hat_ct);
	my $mut_str = sprintf("%-70s %-10s", "Results - Number of Mutator elements:", $mut_ct);
	my $cac_str = sprintf("%-70s %-10s", "Results - Number of CACTA elements:", $cac_ct);
	my $mte_str = sprintf("%-70s %-10s", "Results - Number of MITE elements:", $unc_obj->{mite_count});
	my $unc_str = sprintf("%-70s %-10s", "Results - Number of unclassified TIR elements:", $unc_obj->{unc_count});

	$log->info($tot_str);
	$log->info($tc1_str);
	$log->info($hat_str);
	$log->info($mut_str);
	$log->info($cac_str);
	$log->info($mte_str);
	$log->info($unc_str);

	return (\%gffs, $unc_obj->{mite_index}, $log);
    }
    else {
	say STDERR "\n[WARNING]: No TIR elements were classified. Check input.\n";
    }
    return;
}

sub _classify_tir_families {
    my ($opt, $gffs, $mite_index, $log) = @_;

    my $threads = $opt->{threads} // 1;
    my $hpcov   = $opt->{percentcov} // 50;
    my $hpid    = $opt->{percentid} // 80;
    my $hlen    = $opt->{hitlen} // 80;
    my $debug   = $opt->{debug} // 0;
    my $type    = 'TIR';

    unless ( -d $opt->{outdir} ) {
        make_path( $opt->{outdir}, {verbose => 0, mode => 0771,} );
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

    # update the FASTA IDs and combine them
    $classify_fams_obj->combine_families({ outfiles      => $outfiles, 
                                           annotated_ids => $annot_ids, 
                                           annotated_idx => $mite_index });

    # update the GFF3 IDs and combine them
    $classify_fams_obj->annotate_gff({ annotated_ids => $annot_ids, 
                                       annotated_idx => $mite_index, 
                                       input_gff     => $opt->{ingff}, 
                                       te_type       => 'TIR' });

    
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

S. Evan Staton, C<< <evan at evanstaton.com> >>

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
