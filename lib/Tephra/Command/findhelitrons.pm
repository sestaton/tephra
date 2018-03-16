package Tephra::Command::findhelitrons;
# ABSTRACT: Find Helitons in a genome assembly.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use Cwd           qw(abs_path);
use File::Copy    qw(move);
use File::Basename;
use Tephra -command;
use Tephra::Config::Exe;
use Tephra::Hel::HelSearch;
use Tephra::Classify::Any;
#use Data::Dump::Color;

sub opt_spec {
    return (    
	[ "genome|g=s",           "The genome sequences in FASTA format to search for Helitrons "   ],
	[ "helitronscanner|j=s",  "The HelitronScanner .jar file (configured automatically) "       ],
	[ "gff|o=s",              "The final combined and filtered GFF3 file of Helitrons "         ],
	[ "fasta|f=s",            "The final combined and filtered FASTA file of Helitrons "        ],
	[ "threads|t=i",          "The number of threads to use for BLAST searches (Default: 1)  "  ],
	[ "logfile=s",            "The file to use for logging results in addition to the screen "  ],
	[ "debug",                "Show external command for debugging (Default: no) "              ],
	[ "help|h",               "Display the usage menu and exit. "                               ],
        [ "man|m",                "Display the full manual. "                                       ],
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
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my ($hel_obj, $sf_elem_map) = _run_helitron_search($opt);
    _find_helitron_families($opt, $hel_obj, $sf_elem_map);
}

sub _run_helitron_search {
    my ($opt) = @_;
    
    my $genome   = $opt->{genome};
    my $hscan    = $opt->{helitronscanner};
    my $gff      = $opt->{gff};
    my $debug    = $opt->{debug} // 0;
    my $config   = Tephra::Config::Exe->new->get_config_paths;
    my ($hscanj) = @{$config}{qw(hscanjar)};

    $hscan //= $hscanj;

    #say STDERR "hscandir: $hscan";
    my %opts = (
	genome          => $genome,
	helitronscanner => $hscan,
        gff             => $gff,
        debug           => $debug );

    if (defined $opt->{fasta}) {
	$opts{fasta} = $opt->{fasta};
    }

    my $hel_search = Tephra::Hel::HelSearch->new(%opts);
    my $hel_seqs = $hel_search->find_helitrons;
    my ($sf_elem_map, $hel_obj) = $hel_search->make_hscan_outfiles($hel_seqs);
    
    return ($hel_obj, $sf_elem_map);
}

sub _find_helitron_families {
    my ($opt, $hel_obj, $sf_elem_map) = @_;
 
    my ($name, $path, $suffix) = fileparse($opt->{gff}, qr/\.[^.]*/);
    my $fasta = $opt->{fasta} // File::Spec->catfile($path, $name.'.fasta');
    my $threads = $opt->{threads} // 1;

    my $anno_obj = Tephra::Classify::Any->new(
        fasta   => $hel_obj->{fasta},
        gff     => $opt->{gff},
        threads => $threads,
        outdir  => $path,
    );
 
    my ($logfile, $log);
    if ($opt->{logfile}) {
        $log = $anno_obj->get_tephra_logger($opt->{logfile});
    }
    else {
        my ($gname, $gpath, $gsuffix) = fileparse($opt->{genome}, qr/\.[^.]*/);
        $logfile = File::Spec->catfile( abs_path($gpath), $gname.'_tephra_findhelitrons.log' );
        $log = $anno_obj->get_tephra_logger($logfile);
        say STDERR "\n[WARNING]: '--logfile' option not given so results will be appended to: $logfile.";
    }

    my $blast_report = $anno_obj->process_blast_args;

    if (defined $blast_report) { 
	my $matches = $anno_obj->parse_blast($blast_report);
	my ($fams, $ids, $sfmap, $family_stats) = 
	    $anno_obj->write_families($hel_obj->{fasta}, $matches, $sf_elem_map, 'helitron');
	my $totct = $anno_obj->combine_families($fams, $fasta);
	$anno_obj->annotate_gff($ids, $hel_obj->{gff}, $sf_elem_map);
	
	my ($elemct, $famct, $famtot, $singct) =
	    @{$family_stats}{qw(total_elements families total_in_families singletons)};
	
	$log->info("Results - Number of Helitron families:                         $famct");
	$log->info("Results - Number of Helitron elements in families:             $famtot");
	$log->info("Results - Number of Helitron singleton families/elements:      $singct");
	$log->info("Results - Number of Helitron elements (for debugging):         $elemct");
	$log->info("Results - Number of Helitron elements written (for debugging): $totct");
	
	unlink $_ for keys %$fams;
	unlink @{$hel_obj}{qw(fasta gff)};
    }
    else {
	say "\n[WARNING]: No BLAST hits were found so no Helitron families could be determined.\n";
	move $hel_obj->{fasta}, $fasta or die "\n[ERROR]: move failed: $!\n";
	move $hel_obj->{gff}, $opt->{gff} or die "\n[ERROR]: move failed: $!\n";
    }
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra findhelitrons [-h] [-m]
    -m --man                :   Get the manual entry for a command.
    -h --help               :   Print the command usage.

Required:
    -g|genome               :   The genome sequences in FASTA format to search for Helitrons.. 
    -o|gff                  :   The final combined and filtered GFF3 file of Helitrons.

Options:
    -f|fasta                :   The final combined and filtered FASTA file of Helitrons.
                                (Default name is that same as the GFF3 file except with the ".fasta" extension)
    -d|helitronscanner_dir  :   The HelitronScanner directory containing the ".jar" files and Training Set.
                                This should be configured automatically upon a successful install.
    -t|threads              :   The number of threads to use for BLAST searches (Default: 1).
    --logfile               :   The file to use for logging results in addition to the screen.
    --debug                 :   Show external command for debugging (Default: no).

END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findhelitrons - Find Helitrons in a genome assembly.

=head1 SYNOPSIS    

 tephra findhelitrons -g ref.fas -o ref_helitrons.gff3

=head1 DESCRIPTION

 Find Helitionrs in a reference genome assembly.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for TIR TEs.

=item -o, --gff

 The final combined and filtered GFF3 file of Helitrons.

=back

=head1 OPTIONS

=over 2

=item -f, --fasta

  The final combined and filtered FASTA file of Helitrons.  

=item -d, --helitronscanner_dir

 The HelitronScanner directory. This should not have to be used except by developers as it
 should be configured automatically during the installation.

=item -t, --threads

 The number of threads to use for BLAST searches (Default: 1).

=item --logfile

 The file to use for logging results in addition to the screen.

=item --debug

 Show external command for debugging (Default: no).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
