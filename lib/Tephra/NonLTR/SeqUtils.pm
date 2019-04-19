package Tephra::NonLTR::SeqUtils;

use 5.014;
use Moose;
use File::Find;
use File::Basename;
use File::Path          qw(make_path);
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(system);
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::NonLTR::SeqUtils - Minor sequence utilities for non-LTR finding

=head1 VERSION

Version 0.12.3

=cut

our $VERSION = '0.12.3';
$VERSION = eval $VERSION;

has verbose => ( is => 'ro', isa => 'Bool', predicate  => 'has_verbose', lazy => 1, default => 0 );

sub invert_seq {
    my $self = shift;
    my ($plus_dna_dir, $minus_dna_dir) = @_;

    unless ( -d $minus_dna_dir ) {
        make_path( abs_path($minus_dna_dir), {verbose => 0, mode => 0771,} );
    }

    my @fasfiles;
    find( sub { push @fasfiles, $File::Find::name if -f }, $plus_dna_dir );

    my @revfasfiles;
    for my $file (@fasfiles) {
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);

	open my $in, '<', $file or die "\n[ERROR]: Could not open file: $file";
        my @temp = <$in>;
        close $in;

	shift @temp if $temp[0] =~ /^\>/;
        chomp @temp;
        my $seq = join "", @temp;
        my $revseq = reverse $seq;
        $revseq =~ tr/[A,C,G,T,a,c,g,t]/[T,G,C,A,t,g,c,a]/;
        $revseq =~ s/.{60}\K/\n/g;
        my $outfile = File::Spec->catfile($minus_dna_dir, $name.$suffix);
        open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile";;
	say $out join "\n", ">".$name, $revseq;
	close $out;
	push @revfasfiles, $outfile;
    }

    return \@revfasfiles;
}

sub translate {
    my $self = shift;
    my ($in, $out, $strand) = @_;
    #my $pdir = $self->pdir->absolute->resolve;

    my $name = basename($in);
    #my $config = Tephra::Config::Exe->new->get_config_paths;                                                                               
    #my ($translate) = @{$config}{qw(transcmd)};                                                                                            
    #my $cmd = "$translate -d $in -h $name -p $out";                                                                                        

    my $frame = $strand =~ /forward|plus/i ? 'F' 
	      : $strand =~ /reverse|minus/i ? 'R' 
	      : 0;

    unless (defined $frame) {
	say STDERR "\n[ERROR]: Could not determine frame for translation. Exiting.\n";
	return $frame;
    }

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($transeq) = @{$config}{qw(transeq)};
    #my $cmd = "transeq -frame R -sequence t/test_data/Ha412HOChr01_genome/Ha412HOChr01.fasta -outseq t/test_data/Ha412HOChr01_nonLTRs/b/Ha412HOChr01_rev_trans_trim_clean.faa -auto -trim -clean -reverse"; 
    my $cmd = join q{ }, $transeq, '-frame', $frame, '-sequence', $in, '-outseq', $out, '-trim', '-clean', '-auto';
    say STDERR "CMD: $cmd" if $self->verbose;

    try {
        system($cmd);
        #system([0..5], $transeq, '-frame', $frame, '-sequence', $dna_file, '-outseq', $pep_file, '-trim', '-clean', '-auto');
    }
    catch {
        say STDERR "\n[ERROR]: 'transeq' died. Here is the exception: $_\n";
        exit(1);
    };

    # frame is defined and true, so we return that to say all went well
    return $frame;
}

=head1 AUTHOR

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::NonLTR::SeqUtils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
