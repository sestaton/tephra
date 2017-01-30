package Tephra::Command::findhelitrons;
# ABSTRACT: Find Helitons in a genome assembly.

use 5.014;
use strict;
use warnings;
use File::Basename;
use Tephra -command;
use Tephra::Config::Exe;
use Tephra::Hel::HelSearch;

sub opt_spec {
    return (    
	[ "genome|g=s",           "The genome sequences in FASTA format to search for Helitrons "   ],
	[ "helitronscanner|j=s",  "The HelitronScanner .jar file (configured automatically) "       ],
	[ "outgff|o=s",           "The final combined and filtered GFF3 file of Helitrons "         ],
	[ "outfasta|f=s",         "The final combined and filtered FASTA file of Helitrons "        ],
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
        $self->help;
        exit(0);
    }
    elsif (!$opt->{genome} || !$opt->{outgff}) {
	say STDERR "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $gff = _run_helitron_search($opt);
}

sub _run_helitron_search {
    my ($opt) = @_;
    
    my $genome   = $opt->{genome};
    my $hscan    = $opt->{helitronscanner};
    my $gff      = $opt->{outgff};
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

    if (defined $opt->{outfasta}) {
	$opts{fasta} = $opt->{outfasta};
    }

    my $hel_search = Tephra::Hel::HelSearch->new(%opts);

    my $hel_seqs = $hel_search->find_helitrons;
    $hel_search->make_hscan_outfiles($hel_seqs);
    
    return $gff;
}

sub help {
    print STDERR<<END

USAGE: tephra findhelitrons [-h] [-m]
    -m --man                :   Get the manual entry for a command.
    -h --help               :   Print the command usage.

Required:
    -g|genome               :   The genome sequences in FASTA format to search for Helitrons.. 
    -o|outgff               :   The final combined and filtered GFF3 file of Helitrons.

Options:
    -f|outfasta             :   The final combined and filtered FASTA file of Helitrons.
                                (Default name is that same as the GFF3 file except with the ".fasta" extension)
    -d|helitronscanner_dir  :   The HelitronScanner directory containing the ".jar" files and Training Set.
                                This should be configured automatically upon a successful install.
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

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for TIR TEs.

=item -o, --outgff

 The final combined and filtered GFF3 file of Helitrons.

=back

=head1 OPTIONS

=over 2

=item -f, --outfasta

  The final combined and filtered FASTA file of Helitrons.  

=item -d, --helitronscanner_dir

 The HelitronScanner directory. This should not have to be used except by developers as it
 should be configured automatically during the installation.

=item --debug

  Show external command for debugging (Default: no).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
