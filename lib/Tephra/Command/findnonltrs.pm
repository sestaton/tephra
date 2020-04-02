package Tephra::Command::findnonltrs;
# ABSTRACT: Find non-LTR retrotransposons in a genome assembly.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Cwd           qw(abs_path);
use File::Copy    qw(move);
use File::Path    qw(remove_tree);
use File::Spec;
use File::Basename;
use Tephra -command;
use Tephra::NonLTR::NonLTRSearch;
use Tephra::NonLTR::GFFWriter;
use Tephra::Classify::Any;
use Carp 'croak';
#use Data::Dump::Color;

sub opt_spec {
    return (    
	[ "genome|g=s",    "The genome sequences in FASTA format to search for non-LTR-RTs " ],
	[ "reference|r=s", "The non-masked reference genome for base correction "            ], 
	[ "pdir|d=s",      "The directory to search for HMMs (configured automatically) "    ],
	[ "outdir|d=s",    "The location to place the results "                              ],
	[ "gff|o=s",       "The GFF3 outfile to place the non-LTRs found in <genome> "       ],
	[ "threads|t=i",   "The number of threads to use for BLAST searches (Default: 1)  "  ],
	[ "logfile|l=s",   "The file to use for logging results in addition to the screen "  ],
	[ "verbose|v",     "Display progress for each chromosome (Default: no) "             ],
	[ "help|h",        "Display the usage menu and exit. "                               ],
        [ "man|m",         "Display the full manual. "                                       ],
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
    elsif (!$opt->{genome} || !$opt->{gff}) {
	say STDERR "\n[ERROR]: Required arguments not given.\n";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{genome}) {
        say STDERR "\n[ERROR]: The genome file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
        $self->app->global_options->{help};

    my ($nonltr_obj, $sf_elem_map) = _run_nonltr_search($opt);
    if (defined $nonltr_obj && defined $sf_elem_map) { 
	_find_nonltr_families($opt, $nonltr_obj, $sf_elem_map);
    }
}

sub _run_nonltr_search {
    my ($opt) = @_;
    
    my $genome    = $opt->{genome};
    my $reference = $opt->{reference};
    my $outdir    = $opt->{outdir};
    my $gff       = $opt->{gff};
    my $pdir      = $opt->{pdir} // $ENV{TEPHRA_DIR} // File::Spec->catdir($ENV{HOME}, '.tephra');
    my $verbose   = $opt->{verbose} // 0;
    my $threads   = $opt->{threads} // 1;

    my $nonltr_obj = Tephra::NonLTR::NonLTRSearch->new(
	genome  => $genome,
	outdir  => $outdir,
	pdir    => $pdir,
	threads => $threads,
	verbose => $verbose );
    
    my ($genomedir, $outputdir) = $nonltr_obj->find_nonltrs;

    if (defined $genomedir && defined $outputdir) {
	$genome = defined $reference ? $reference : $genome;
	my $gff_obj = Tephra::NonLTR::GFFWriter->new(
	    genome   => $genome,
	    fastadir => $genomedir,
	    outdir   => $outputdir,
	    gff      => $gff );
	
	my ($obj, $sf_elem_map) = $gff_obj->write_gff;
	
	return ($obj, $sf_elem_map);
    }
    else {
	say STDERR "\n[WARNING]: No non-LTR elements were found so none will be reported.\n";
	return (undef, undef);
    }
}

sub _find_nonltr_families {
    my ($opt, $obj, $sf_elem_map) = @_;

    my ($name, $path, $suffix) = fileparse($opt->{gff}, qr/\.[^.]*/);
    my $fasta = $opt->{fasta} // File::Spec->catfile($path, $name.'.fasta');
    my $threads = $opt->{threads} // 1;

    my $anno_obj = Tephra::Classify::Any->new(
	gff     => $opt->{gff},
	threads => $threads,
	type    => 'non-LTR',
     );
    
    my $seqstore = $anno_obj->store_seq($obj->{fasta});

    my ($logfile, $log);
    if ($opt->{logfile}) {
        $log = $anno_obj->get_tephra_logger($opt->{logfile});
    }
    else {
        my ($gname, $gpath, $gsuffix) = fileparse($opt->{genome}, qr/\.[^.]*/);
        $logfile = File::Spec->catfile( abs_path($gpath), $gname.'_tephra_findnonltrs.log' );
        $log = $anno_obj->get_tephra_logger($logfile);
        say STDERR "\n[WARNING]: '--logfile' option not given so results will be appended to: $logfile.\n";
    }

    my (%families, %family_files, %family_map);
    for my $elem (keys %$sf_elem_map) {
        push @{$families{ $sf_elem_map->{$elem} }}, $elem;
    }
 
    for my $fam (keys %families) { 
        my $famfile = File::Spec->catfile( abs_path($path), $fam.'_elements.fasta' );
        $family_files{$fam} = $famfile;

        open my $out, '>', $famfile or die "\n[ERROR]: Could not open file: $famfile\n";
        
        for my $elem (@{$families{$fam}}) {
            if (exists $seqstore->{$elem}) {
                say $out join "\n", ">$elem", $seqstore->{$elem};
            }
            else {
                croak "\n[ERROR]: $elem not found in store. Exiting.";
            }
        }
	close $out;
    }
    #say "===> FAMILY_FILES";
    #dd \%family_files;#

    my ($elemct, $famct, $singct, $famtot) = (0, 0, 0, 0);
    for my $family (sort keys %family_files) {
	my $blast_report = $anno_obj->process_blast_args($family_files{$family});

	my $matches;
	if (defined $blast_report) {
	    $matches = $anno_obj->parse_blast($blast_report);
	}
	else {
	    $matches = {};
	}
	
	my ($fams, $ids, $family_stats) = 
	    $anno_obj->write_families($family_files{$family}, $matches, $sf_elem_map, $famct, $singct, $famtot);

	#say join q{ }, "DEBUG:", "FAMILY", "TOTAL_ELEMENTS", "FAMILIES", "TOTAL_IN_FAMILIES", "SINGLETONS";
	#say join q{ }, "         $family", $family_stats->{total_elements}, $family_stats->{families}, 
	#$family_stats->{total_in_families}, $family_stats->{singletons};

	$elemct += $family_stats->{total_elements};
	$famct  += $family_stats->{families};
	$famtot += $family_stats->{total_in_families};
	$singct += $family_stats->{singletons};

	# add files to objects
	$family_map{$family} = { IDS => $ids, FAMS => $fams };

	#else {
	    #say STDERR "\n[WARNING]: No BLAST hits were found so no non-LTR families could be determined.\n";
            #move $obj->{fasta}, $fasta or die "\n[ERROR]: move failed: $!\n";
            #move $obj->{gff}, $opt->{gff} or die "\n[ERROR]: move failed: $!\n";  
	#}
    }

    my $totct = $anno_obj->combine_families(\%family_map, $fasta);
    $anno_obj->annotate_gff(\%family_map, $obj->{gff}, $sf_elem_map);
    
    $log->info("Results - Number of non-LTR families:                         $famct");
    $log->info("Results - Number of non-LTR elements in families:             $famtot");
    $log->info("Results - Number of non-LTR singleton families/elements:      $singct");
    $log->info("Results - Number of non-LTR elements (for debugging):         $elemct");
    $log->info("Results - Number of non-LTR elements written (for debugging): $totct");
    
    unlink $_ for values %family_files;
    unlink @{$obj}{qw(fasta gff)};
    for my $fam (keys %family_map) {
	for my $fam_file (keys %{$family_map{$fam}{FAMS}}) {
	    unlink $fam_file;
	}
    }

    # clean up analysis directories after searching for families
    for my $dir (@{$obj->{fasta_dirs}}) {
	remove_tree( $dir, { safe => 1 } );
    }

    return;
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra findnonltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for non-LTR-RTs. 
    -o|gff        :   The GFF3 outfile to place the non-LTRs found in <genome>.

Options:
    -r|reference  :   The non-masked reference genome for base correction.
    -d|outdir     :   The location to place the results.
    -p|pdir       :   Location of the HMM models (Default: configured automatically).
    -t|threads    :   The number of threads to use for BLAST searches (Default: 1).
    -v|verbose    :   Display progress for each chromosome (Default: no).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findnonltrs - Find non-LTRs retrotransposons in a genome assembly.

=head1 SYNOPSIS    

 tephra findnonltrs -g genome.fas -d ref_nonltrs_results -o genome_nonltrs.gff3

=head1 DESCRIPTION
 
 Find non-LTR retrotransposons in a reference genome, classify them into known superfamilies, 
 and generate a GFF file showing their locations and properties.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

The genome sequences in FASTA format to search for non-LTR-RTs.

=item o, --gff

 The GFF3 outfile to place the non-LTRs found in <genome>.

=back

=head1 OPTIONS

=over 2

=item -r, --reference

 The non-masked reference genome for base correction. If a masked genome was used to search for non-LTRs,
 it is desirable to correct the bases using the original non-masked reference sequence.

=item -o, --outdir

 The directory name to place the resulting GFF file (and run the analysis).

=item -d, --pdir

 The directory to search for HMMs. This should be configured automatically during installation and this option should only have to be used by developers.

=item -t, --threads

 The number of threads to use for BLAST searches (Default: 1).

=item -v, --verbose

 Display progress for each chromosome (Default: no).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
