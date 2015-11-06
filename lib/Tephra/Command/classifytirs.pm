package Tephra::Command::classifytirs;
# ABSTRACT: Classify TIR transposons into superfamilies.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::Classify::TIRSfams;
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use File::Basename;
use File::Spec;

sub opt_spec {
    return (    
	[ "genome|g=s",   "The genome sequences in FASTA format to search for TIRs "   ],
	[ "gff|f=s",      "The GFF3 file of TIR TEs in <genome> "                      ],
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

    my $some = _classify_tir_predictions($opt);
}

sub _classify_tir_predictions {
    my ($opt) = @_;

    my $classify_obj = Tephra::Classify::TIRSfams->new( 
	genome   => $opt->{genome}, 
	gff      => $opt->{gff} 
    );

    my ($header, $features) = $classify_obj->collect_gff_features($opt->{gff});

    my $all_ct = (keys %$features);
    $classify_obj->find_tc1_mariner($features, $header);
    my $tc1_ct = (keys %$features);
    $classify_obj->find_hat($features, $header);
    my $hat_ct = (keys %$features);
    $classify_obj->find_mutator($features, $header);
    my $mut_ct = (keys %$features);
    $classify_obj->find_cacta($features, $header);
    my $cacta_ct = (keys %$features);
    $classify_obj->write_unclassified_tirs($features, $header);
    my $rem_ct = (keys %$features);

    say STDERR join "\t", "all", "after_tc1", "after_hat", "after_mut", "after_cacta", "after_rem";
    say STDERR join "\t", $all_ct, $tc1_ct, $hat_ct, $mut_ct, $cacta_ct, $rem_ct;
}

sub help {
    print STDERR<<END

USAGE: tephra classifytirs [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for TIR TEs. 
    -f|gff        :   The GFF3 file of LTR-RTs in <--genome>.

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra classifytirs - Classify TIR transposons into superfamilies.

=head1 SYNOPSIS    

 tephra findltrs -g ref.fas -f ref_tephra_tirs.gff3

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

=back

=head1 OPTIONS

=over 2

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
