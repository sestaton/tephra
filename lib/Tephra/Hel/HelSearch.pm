package Tephra::Hel::HelSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Sort::Naturally;
use IPC::System::Simple qw(system EXIT_ANY);
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;
use namespace::autoclean;

with 'Tephra::Role::Run::HelitronScanner';

=head1 NAME

Tephra::Hel::HelSearch - Find Helitrons in a reference genome

=head1 VERSION

Version 0.03.5

=cut

our $VERSION = '0.03.5';
$VERSION = eval $VERSION;

sub find_helitrons {
    my $self = shift;
    
    my $genome = $self->genome->absolute;
    my $jar    = $self->helitronscanner;
    my $gff    = $self->outfile;

    my (%scanh_cmd, %scant_cmd, %pair_cmd, %draw_cmd);
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
 
    my $g_headlcvs = File::Spec->catfile($path, $name.'_hscan_head.lcvs');
    my $g_taillcvs = File::Spec->catfile($path, $name.'_hscan_tail.lcvs');
    my $g_paired   = File::Spec->catfile($path, $name.'_hscan_paired.txt');
    my $g_helname  = File::Spec->catfile($path, $name.'_tephra_hscan_helitrons');
    my $full_hels  = $g_helname.'.hel.fa';

    #my $jar  = File::Spec->catfile($hscan_dir, "HelitronScanner.jar");
    my $parent = $jar->parent->parent; # unfortunately, the dist does not unpack in a separate dir
    my $lcvs = File::Spec->catfile($parent, 'TrainingSet', 'head.lcvs');
    my $rcvs = File::Spec->catfile($parent, 'TrainingSet', 'tail.lcvs');

    my @scanh_opts = qw(-g -lf -o -tl -buffer_size);
    my @scanh_args = ($genome, $lcvs, $g_headlcvs, '10', '1000000');
    @scanh_cmd{@scanh_opts} = @scanh_args;

    my @scant_opts = qw(-g -lf -o -tl -buffer_size);
    my @scant_args = ($genome, $rcvs, $g_taillcvs, '10', '1000000');
    @scant_cmd{@scant_opts} = @scant_args;

    my @pair_opts = qw(-hs -ts -o);
    my @pair_args = ($g_headlcvs, $g_taillcvs, $g_paired);
    @pair_cmd{@pair_opts} = @pair_args;

    my @draw_opts = qw(-p -g -o -ext5 -ext3);
    my @draw_args = ($g_paired, $genome, $g_helname, '100', '100');
    @draw_cmd{@draw_opts} = @draw_args;
    
    $self->run_hscan_headtail(\%scanh_cmd, $jar, 'scanHead');
    $self->run_hscan_headtail(\%scant_cmd, $jar, 'scanTail');
    $self->run_hscan_pair(\%pair_cmd, $jar);
    $self->run_hscan_draw(\%draw_cmd, $jar);

    return $full_hels;
}

sub make_hscan_gff {
    my $self = shift;
    my ($helitrons) = @_;
    my $gff = $self->outfile;
    my $genome = $self->genome;
    open my $out, '>', $gff or die "\nERROR: Could not open file: $gff\n";
    
    my %refs;
    my $gkseq = Bio::DB::HTS::Kseq->new($genome);
    my $giter = $gkseq->iterator;

    while (my $gseqs = $giter->next_seq) {
	my $name = $gseqs->name;
	my $seq  = $gseqs->seq;
	$refs{$name} = length($seq);
    }

    my $header = "##gff-version 3";
    say $out $header;
    for my $ref (nsort keys %refs) {
	say $out join q{ }, "##sequence-region", $ref, '1', $refs{$ref};
    }
    
    my %strand = ( forward => '+', reverse => '-' );
    
    my ($id, $seq, %hel);
    my $helct = 0;
    open my $hin, '<', $helitrons or die "\nERROR: Could not open file: $helitrons\n";
    while (($id, $seq) = $self->read_seq(\*$hin)) {
	$helct++;
	my ($ref, $start, $stop) = ($id =~ /(^\S+)_\#SUB_(\d+)-(\d+)/);
	my ($str) = ($id =~ /\[(forward|reverse)\]/);
	my $strand = $strand{$str};

	# seqid source type start end score strand phase attribs
	my $gff_str;
	if ($start > $stop && $strand eq '-') {
	    $gff_str = join "||", $ref, 'HelitronScanner', 'helitron', $stop, $start, '.', 
	        $strand, '.', "ID=helitron$helct;Ontology_term=SO:0000544";
	}
	else {
	    $gff_str = join "||", $ref, 'HelitronScanner', 'helitron', $start, $stop, '.',
                $strand, '.', "ID=helitron$helct;Ontology_term=SO:0000544";
	}
	push @{$hel{$ref}}, $gff_str;
    }
    close $hin;

    for my $ref (nsort keys %hel) {
	for my $feature (@{$hel{$ref}}) {
	    my @feats = split /\|\|/, $feature;
	    say $out join "\t", @feats;
	}
    }
    close $out;
}

sub read_seq {
    my $self = shift;
    my ($fh) = @_;
    
    local $/ = "\n>";
    return unless my $entry = $fh->getline;
    chomp $entry;

    my ($id, $seq) = split /\n/, $entry, 2;
    defined $id && $id =~ s/>//g;
    defined $seq && $seq =~ s/>//g;

    return ($id, $seq);
}

=head1 AUTHOR

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Hel::HelSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
