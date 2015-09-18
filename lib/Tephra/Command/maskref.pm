package Tephra::Command::maskref;
# ABSTRACT: Mask a reference genome with transposons.

use 5.010;
use strict;
use warnings;
use Tephra -command;
use Tephra::MaskRef;
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(system);
use Capture::Tiny       qw(:all);
use File::Basename;
use File::Spec;

sub opt_spec {
    return (    
	[ "genome|g=s",   "The genome sequences in FASTA format to search for LTR-RTs " ],
	[ "repeatdb|d=s", "The database of repeat sequences to use for masking "        ],
	[ "clean",        "Clean up the index files (Default: yes) "                    ],
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
    elsif (!$opt->{genome} || !$opt->{repeatdb}) {
	say "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    exit(0) if $self->app->global_options->{man} ||
	$self->app->global_options->{help};

    my $gff = _run_masking($opt);
}

sub _run_masking {
    my ($opt) = @_;
    
    my $genome   = $opt->{genome};
    my $repeatdb = $opt->{repeatdb};
    my $clean    = defined $opt->{clean} ? $opt->{clean} : 0;

    my $mask_obj = Tephra::MaskRef->new( 
	genome   => $genome, 
	repeatdb => $repeatdb, 
	clean    => $clean 
    );
    
    my $masked_ref = $mask_obj->mask_reference;
}

sub help {
    print STDERR<<END

USAGE: tephra maskref [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -g|genome     :   The genome sequences in FASTA format to search for TIR TEs. 
    -d|repeatdb   :   The database of repeat sequences to use for masking.

Options:
    -c|clean      :   Clean up the index files (Default: yes).

END
}


1;
__END__

=pod

=head1 NAME
                                                                       
 tephra maskref - Mask a reference genome with transposons.

=head1 SYNOPSIS    

 tephra maskref -g ref.fas -d repeatdb.fas

=head1 DESCRIPTION

 Mask a reference genome with one type of transposons to reduce false positives, and
 search time, in subsequent searches for other transposon types.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -g, --genome

 The genome sequences in FASTA format to search for TIR TEs.

=item -d, --repeatdb

 The database of repeat sequences to use for masking.

=back

=head1 OPTIONS

=over 2

=item -c, --clean

 Clean up the index files (Default: yes).

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
