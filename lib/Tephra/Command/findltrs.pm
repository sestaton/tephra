package Tephra::Command::findltrs;
# ABSTRACT: Find LTR retrotransposons in a genome assembly.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Capture::Tiny qw(capture_merged);
use File::Path    qw(remove_tree);
use File::Copy    qw(copy);
use Cwd           qw(getcwd);
use File::Find;
use File::Basename;
use Tephra -command;
use Tephra::Config::Exe;
use Tephra::Config::Reader;
use Tephra::LTR::LTRSearch;
use Tephra::LTR::LTRRefine;
#use Data::Dump::Color;

sub opt_spec {
    return (    
	[ "config|c=s",  "The Tephra LTR option configuration file "                      ],
	[ "outfile|o=s", "The final combined and filtered GFF3 file of LTR-RTs "          ],
	[ "logfile|l=s", "The file to use for logging results in addition to the screen " ],
	[ "index|i=s",   "The suffixerator index to use for the LTR search "              ],
	[ "threads|t=i", "The number of threads to use for the LTR search (Default: 1) "  ],
	[ "help|h",      "Display the usage menu and exit. "                              ],
        [ "man|m",       "Display the full manual. "                                      ],
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
    elsif (!$opt->{config}) {
	say STDERR "\nERROR: Required arguments not given.\n";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{config}) { 
	say STDERR "\nERROR: '--config' file given but does not appear to exist. Check input.\n";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my ($global_opts, $search_config, $relaxed_gff, $strict_gff) = _run_ltr_search($opt);
    my $some = _refine_ltr_predictions($global_opts, $search_config, $relaxed_gff, $strict_gff, $opt);
}

sub _refine_ltr_predictions {
    my ($global_opts, $search_config, $relaxed_gff, $strict_gff, $opt) = @_;

    my %refine_opts = (
	genome => $global_opts->{genome}, 
    );

    $refine_opts{remove_dup_domains} = $search_config->{findltrs}{dedup} =~ /yes/i ? 1 : 0;
    $refine_opts{remove_tnp_domains} = $search_config->{findltrs}{tnpfilter} =~ /yes/i ? 1 : 0;
    $refine_opts{domains_required}   = $search_config->{findltrs}{domains_required} =~ /yes/i ? 1 : 0;
    $refine_opts{outfile} = $opt->{outfile} if $opt->{outfile};
    $refine_opts{logfile} = $opt->{logfile} if $opt->{logfile};
    
    my $refine_obj = Tephra::LTR::LTRRefine->new(%refine_opts);

    if ($relaxed_gff && $strict_gff) {
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
	
	unlink $relaxed_gff, $strict_gff;
    }
    elsif ($relaxed_gff && !$strict_gff) {
	say STDERR "\nWARNING: No LTR retrotransposons were found under strict conditions. ".                  
            "Skipping refinement step.\n";
	$refine_obj->sort_features({ gff               => $relaxed_gff,
                                     combined_features => undef });
	
	unlink $relaxed_gff;
    }
    else {
	say STDERR "\nWARNING: No LTR retrotransposons were found with the given parameters.\n";
    }
}
    
sub _run_ltr_search {
    my ($opt) = @_;
    
    my @indexfiles;
    if (defined $opt->{index}) {
	my ($name, $path, $suffix) = fileparse($opt->{index}, qr/\.[^.]*/);
	my @files;
	for my $suf ('.des', '.lcp', '.llv', '.md5', '.prj', '.sds', '.suf')  {
	    push @files, $opt->{index}.$suf;
	}
	
	my $matchstr = join "|", @files;
	find( sub { push @indexfiles, $File::Find::name if -f and /$matchstr/ }, $path );
    }
    
    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($tephra_hmmdb, $tephra_trnadb) = @{$config}{qw(hmmdb trnadb)};

    my $config_obj    = Tephra::Config::Reader->new( config => $opt->{config} );
    my $search_config = $config_obj->get_configuration;
    my $global_opts   = $config_obj->get_all_opts($search_config);
   
    my $using_tephra_db = 0;
    if ($global_opts->{hmmdb} =~ /TephraDB/) { 
	$using_tephra_db = 1;
	my $dir = getcwd();
	my $tmpiname  = 'tephra_transposons_hmmdb_XXXX';
	my $tmp_hmmdbfh = File::Temp->new( TEMPLATE => $tmpiname,
					   DIR      => $dir,
					   SUFFIX   => '.hmm',
					   UNLINK   => 0);
	my $tmp_hmmdb = $tmp_hmmdbfh->filename;
	copy $tephra_hmmdb, $tmp_hmmdb or die "Copy failed: $!";
	$global_opts->{hmmdb} = $tmp_hmmdb;
    }

    my $ltr_search_obj = Tephra::LTR::LTRSearch->new( 
	genome  => $global_opts->{genome},
	hmmdb   => $global_opts->{hmmdb},
	trnadb  => $global_opts->{trnadb}, 
	clean   => $global_opts->{clean},
	debug   => $global_opts->{debug},
	threads => $global_opts->{threads},
	logfile => $global_opts->{logfile},
    );

    unless (defined $opt->{index} && @indexfiles == 7) {
	my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
	$opt->{index} = $global_opts->{genome}.'.index';
    
	my @suff_args = qq(-db $global_opts->{genome} -indexname $opt->{index} -tis -suf -lcp -ssp -sds -des -dna);
	$ltr_search_obj->create_index(\@suff_args, $global_opts->{logfile});
    }

    #my ($strict_gff, $strict_models) = 
    my $strict_gff =
	$ltr_search_obj->ltr_search({ config => $search_config, index => $opt->{index}, mode => 'strict'  });
    #my ($relaxed_gff, $relaxed_models) = 
    my $relaxed_gff =
	$ltr_search_obj->ltr_search({ config => $search_config, index => $opt->{index}, mode => 'relaxed' });
    #remove_tree($strict_models, { safe => 1 });
    #remove_tree($relaxed_models, { safe => 1 });
    unlink $global_opts->{hmmdb} if $using_tephra_db; # this is just a temp file to keep ltrdigest from crashing

    return ($global_opts, $search_config, $relaxed_gff, $strict_gff);
}

sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra findltrs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -c|config     :   The Tephra configuration file.

Options:
    -o|outfile    :   The final combined and filtered GFF3 file of LTR-RTs.
    -i|index      :   The suffixerator index to use for the LTR search.
    -t|threads    :   The number of threads to use for the LTR search (Default: 1).

END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra findltrs - Find LTR retrotransposons in a genome assembly.

=head1 SYNOPSIS    

 tephra findltrs --config tephra_ltrs_conf.yml

=head1 DESCRIPTION
 
 Search a reference genome and find LTR-RTs under relaxed and strict conditions (more on
 this later...), filter all predictions by a number of criteria, and generate a consensus
 set of the best scoring elements.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -c, --config

 The Tephra configuration file.

=back

=head1 OPTIONS

=over 2

=item -o, --outfile

 The final combined and filtered GFF3 file of LTR-RTs.

=item -i, --index

 The suffixerator index to use for the LTR search. 

=item -t, --threads

 The number of threads to use for the LTR search. Specifically, the number
 of processors to use by the 'gt' program (Default: 1).

=item -t, --trnadb

 The file of tRNA sequences in FASTA format to search for PBS.

=item -d, --hmmdb

 The HMM db in HMMERv3 format to search for coding domains.

=item --mintsd

 The minimum TSD length (Default: 4).

=item --maxtsd

 The maximum TSD length (Default: 6).

=item --minlenltr

 The minimum LTR length (Default: 100).

=item --maxlenltr

  The maximum LTR length (Default: 6000).

=item --mindistltr

 The minimum LTR element length (Default: 1500).

=item --maxdistltr

 The maximum LTR element length (Default: 25000).

=item --overlaps 

 Keep 'all', 'best', or 'no' overlapping LTR-RT predictions (Default: best).

=item -e, -- pdomevalue

 Protein domain match threshold for pHMM matches with HMMER (Default: 10E-6).

=item -m, --pdomcutoff

 Protein domain match cutoff method for pHMM matches with HMMER. Options are, 'GA' or gathering
 method, 'TC' or trusted cutoff, or 'NONE' (Default: NONE).

=item -r, --dedup

 Discard elements with duplicate coding domains (Default: no).

=item --tnpfilter

 Discard elements containing transposase domains (Default: no).

=item --clean

 Clean up the index files (Default: yes).

=item --debug

 Show external commands for debugging (Default: no).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
