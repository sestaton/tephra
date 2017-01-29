package Tephra::Command::classifytirs;
# ABSTRACT: Classify TIR transposons into superfamilies.

use 5.014;
use strict;
use warnings;
use Tephra -command;
use Tephra::Classify::TIRSfams;

sub opt_spec {
    return (    
	[ "genome|g=s",   "The genome sequences in FASTA format to search for TIRs "   ],
	[ "gff|f=s",      "The GFF3 file of TIR TEs in <genome> "                      ],
	[ "outfile|o=s",  "The final combined and filtered GFF3 file of TIRs "         ],
	[ "help|h",       "Display the usage menu and exit. "                          ],
        [ "man|m",        "Display the full manual. "                                  ],
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
    elsif (!$opt->{genome} || !$opt->{gff} || !$opt->{outfile}) {
	say STDERR "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{genome} || ! -e $opt->{gff}) {
	say STDERR "\nERROR: One or more of the required files does not exist.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $some = _classify_tir_predictions($opt);
}

sub _classify_tir_predictions {
    my ($opt) = @_;

    my $classify_obj = Tephra::Classify::TIRSfams->new( 
	genome   => $opt->{genome}, 
	gff      => $opt->{gff},
	outfile  => $opt->{outfile}
    );

    my $index = $opt->{genome}.'.fai';
    unless (-e $index) {
	$classify_obj->index_ref;
    }
    my ($header, $features) = $classify_obj->collect_gff_features($opt->{gff});

    #my (@fastas, @gffs);
    my $all_ct = (keys %$features);
    my ($tcmoutfile, $tcmfas) = $classify_obj->find_tc1_mariner($features, $header);
    say join q{ }, $tcmoutfile, $tcmfas;
    my $tc1_ct = (keys %$features);
    my ($hatoutfile, $hatfas) = $classify_obj->find_hat($features, $header);
    my $hat_ct = (keys %$features);
    my ($mutoutfile, $mutfas) = $classify_obj->find_mutator($features, $header);
    my $mut_ct = (keys %$features);
    my ($cacoutfile, $cacfas) = $classify_obj->find_cacta($features, $header);
    my $cacta_ct = (keys %$features);
    my ($uncoutfile, $uncfas) = $classify_obj->write_unclassified_tirs($features, $header);
    my $rem_ct = (keys %$features);

    my @fastas = grep { defined && /\.fasta$/ } ($tcmfas, $hatfas, $mutfas, $cacfas, $uncfas);
    my @gffs = grep { defined && /\.gff.*$/ } ($tcmoutfile, $hatoutfile, $mutoutfile, $cacoutfile, $uncoutfile);

    if (@fastas && @gffs) {
	my %outfiles = (
	    fastas => \@fastas,
	    gffs   => \@gffs
	 );

	$classify_obj->write_combined_output(\%outfiles);
	
	say STDERR join "\t", "all", "after_tc1", "after_hat", "after_mut", "after_cacta", "after_rem";
	say STDERR join "\t", $all_ct, $tc1_ct, $hat_ct, $mut_ct, $cacta_ct, $rem_ct;
    }
    else {
	say STDERR "\nWARNING: No TIR elements were classified. Check input.\n";
    }
}
    
sub help {
    print STDERR<<END

USAGE: tephra classifytirs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for TIR TEs. 
    -f|gff        :   The GFF3 file of LTR-RTs in <--genome>.
    -o|outfile    :   The final combined and filtered GFF3 file of TIRs.

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra classifytirs - Classify TIR transposons into superfamilies.

=head1 SYNOPSIS    

 tephra findltrs -g ref.fas -f ref_tephra_tirs.gff3 -o ref_tephra_tirs_classified.gff3

=head1 DESCRIPTION
 
 This subcommand takes a GFF3 as input from Tephra and classifies the TIR TEs into
 superfamilies.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -f, --gff

 The GFF3 file of LTR-RTs in <--genome> as output by the 'tephra findtirs' command.

=item -o, --outfile

 The final combined and filtered GFF3 file of TIRs.

=back

=head1 OPTIONS

=over 2

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
