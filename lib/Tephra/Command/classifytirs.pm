package Tephra::Command::classifytirs;
# ABSTRACT: Classify TIR transposons into superfamilies.

use 5.014;
use strict;
use warnings;
use Cwd qw(abs_path);
use File::Basename;
use Tephra -command;
use Tephra::Classify::TIRSfams;
#use Log::Any qw($log);
#use Data::Dump::Color;

sub opt_spec {
    return (    
	[ "genome|g=s",   "The genome sequences in FASTA format to search for TIRs "       ],
	[ "gff|i=s",      "The GFF3 file of TIR TEs in <genome> "                          ],
	[ "outfile|o=s",  "The final combined and filtered GFF3 file of TIRs "             ],
	[ "logfile=s",    "The file to use for logging results in addition to the screen " ],
	[ "help|h",       "Display the usage menu and exit. "                              ],
        [ "man|m",        "Display the full manual. "                                      ],
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

    my %classify_opts = (
	genome   => $opt->{genome},
        gff      => $opt->{gff},
        outfile  => $opt->{outfile}
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
    my ($header, $features) = $classify_obj->collect_gff_features($opt->{gff});

    my $all_ct = (keys %$features);
    my ($tcmoutfile, $tcmfas, $tc1_ct) = $classify_obj->find_tc1_mariner($features, $header, $index, $log);
    my ($hatoutfile, $hatfas, $hat_ct) = $classify_obj->find_hat($features, $header, $index, $log);
    my ($mutoutfile, $mutfas, $mut_ct) = $classify_obj->find_mutator($features, $header, $index, $log);
    my ($cacoutfile, $cacfas, $cac_ct) = $classify_obj->find_cacta($features, $header, $index, $log);
    my ($uncoutfile, $uncfas, $unc_ct) = $classify_obj->write_unclassified_tirs($features, $header, $index, $log);

    my @fastas = grep { defined && /\.fasta$/ } ($tcmfas, $hatfas, $mutfas, $cacfas, $uncfas);
    my @gffs   = grep { defined && /\.gff3$/  } ($tcmoutfile, $hatoutfile, $mutoutfile, $cacoutfile, $uncoutfile);

    if (@fastas && @gffs) {
	my %outfiles = (
	    fastas => \@fastas,
	    gffs   => \@gffs
	 );

	$classify_obj->write_combined_output(\%outfiles);
	
	$log->info("Results - Total number of TIR elements:                   $all_ct");
	$log->info("Results - Number of Tc1-Mariner elements:                 $tc1_ct");
	$log->info("Results - Number of hAT elements:                         $hat_ct");
	$log->info("Results - Number of Mutator elements:                     $mut_ct");
	$log->info("Results - Number of CACTA elements:                       $cac_ct");
	$log->info("Results - Number of unclassified TIR elements:            $unc_ct");
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
    -i|gff        :   The GFF3 file of TIRs in <--genome>.
    -o|outfile    :   The final combined and filtered GFF3 file of TIRs.

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra classifytirs - Classify TIR transposons into superfamilies.

=head1 SYNOPSIS    

 tephra findltrs -g ref.fas -i ref_tephra_tirs.gff3 -o ref_tephra_tirs_classified.gff3

=head1 DESCRIPTION
 
 This subcommand takes a GFF3 as input from Tephra and classifies the TIR TEs into
 superfamilies.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format used to search for LTR-RTs.

=item -i, --gff

 The GFF3 file of TIRs in <--genome> as output by the 'tephra findtirs' command.

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
