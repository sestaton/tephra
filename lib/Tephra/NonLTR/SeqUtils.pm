package Tephra::NonLTR::SeqUtils;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Temp          qw(tempfile);
use File::Copy          qw(copy move);
use File::Path          qw(make_path);
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(EXIT_ANY system);
use Log::Any            qw($log);
use Bio::DB::HTS::Kseq;
use Tephra::Config::Exe;
use namespace::autoclean;

=head1 NAME

Tephra::NonLTR::SeqUtils - Minor sequence utilities for non-LTR finding

=head1 VERSION

Version 0.13.1

=cut

our $VERSION = '0.13.1';
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
    my ($genome_dir, $in, $out, $strand) = @_;

    my $tmp_file = $in.'.tmp';
    my $frame = $strand =~ /forward|plus/i ? 'F' 
	      : $strand =~ /reverse|minus/i ? 'R' 
	      : 0;
    $frame = 'F';

    unless (defined $frame) {
	say STDERR "\n[ERROR]: Could not determine frame for translation. Exiting.\n";
	return $frame;
    }

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($transeq) = @{$config}{qw(transeq)};
    my @cmd = ($transeq, '-frame', $frame, '-sequence', $in, '-outseq', $tmp_file, '-trim', 'yes', '-auto');
    say STDERR 'CMD: ', join q{ }, @cmd if $self->verbose;

    my $EXIT = system(EXIT_ANY, @cmd);

    if ($EXIT) { # non-zero
        $log->error("'transeq' failed with exit value: $EXIT. Here is the output: \n");
        exit(1);
    }

    ##TODO: Add a check here if the output is defined/exists, and return result appropriately
    my $kseq = Bio::DB::HTS::Kseq->new($tmp_file);
    my $iter = $kseq->iterator;

    say STDERR "=====> Writing translated output: $out" if $self->verbose;
    open my $outfh, '>', $out or die "\n[ERROR]: Could not open file: $out\n";
    my $filename = basename($in);

    while (my $seqio = $iter->next_seq) {
	my $id = $seqio->name;
	my $seq = $seqio->seq;
	$seq =~ s/.{60}\K/\n/g;
	my ($frame) = ($id =~ /(_\d+$)/);
	my $newid = $filename.$frame;
	say $outfh join "\n", ">".$newid, $seq;
    }
    close $outfh;
    unlink $tmp_file;

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
