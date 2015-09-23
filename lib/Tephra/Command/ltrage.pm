package Tephra::Command::ltrage;
# ABSTRACT: Calculate the age distribution of LTR retrotransposons.

use 5.010;
use strict;
use warnings;
use File::Path qw(make_path);
use Tephra -command;
use Tephra::LTR::LTRStats;

sub opt_spec {
    return (
	[ "genome|g=s",   "The genome sequences in FASTA format used to search for LTR-RTs "   ],
	[ "gff|f=s",      "The GFF3 file of LTR-RTs in <genome> "    ],
	[ "outdir|o=s",   "The output directory for placing categorized elements " ],
	[ "threads|t=i",  "The number of threads to use for clustering coding domains " ],
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
    elsif (!$opt->{genome} || !$opt->{gff}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
}

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $some = _calculate_ltr_stats($opt);
}

sub _calculate_ltr_stats {
    my ($opt) = @_;

    my $genome   = $opt->{genome};
    my $gff      = $opt->{gff};
    my $outdir   = $opt->{outdir};
    my $threads  = defined $opt->{threads} ? $opt->{threads} : 1;

    unless ( -d $outdir ) {
	make_path( $outdir, {verbose => 0, mode => 0771,} );
    }

    my $stats_obj = Tephra::LTR::LTRStats->new(
	genome   => $genome,
	gff      => $gff,
	outdir   => $outdir,
    );

    my ($header, $features) = $stats_obj->collect_ltr_features;

}

sub help {
    print STDERR<<END

  USAGE: tephra ltrage [-h] [-m]
      -m --man      :   Get the manual entry for a command.
      -h --help     :   Print the command usage.
    
  Required:
      -g|genome     :   The genome sequences in FASTA format used to search for LTR-RTs.
      -d|repeatdb   :   The file of repeat sequences in FASTA format to use for classification.
      -f|gff        :   The GFF3 file of LTR-RTs in <--genome>.
      -o|outdir     :   The output directory for placing categorized elements.
    
  Options:
      -t|threads    :   The number of threads to use for clustering coding domains (Default: 1).
     
END
}


1;
__END__

=pod
=head1 NAME
                                                                       
 tephra classifyltrs - Classify LTR retrotransposons into superfamilies 
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
=item -h, --help
 Print a usage statement. 
=item -m, --man
 Print the full documentation.
=back
=cut
