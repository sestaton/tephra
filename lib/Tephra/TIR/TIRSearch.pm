package Tephra::TIR::TIRSearch;

use 5.010;
use Moose;
use Cwd;
use Bio::SeqIO;
use File::Spec;
use File::Find;
use File::Copy;
use File::Basename;
use IPC::System::Simple qw(system EXIT_ANY);
use Sort::Naturally;
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;
use namespace::autoclean;

with 'Tephra::Role::Run::GT',
     'Tephra::Role::GFF';

=head1 NAME

Tephra::TIR::TIRSearch - Find TIR transposons in a reference genome

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

sub tir_search {
    my $self = shift;
    my ($index) = @_;
    
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my (%suf_args, %tirv_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $gff = File::Spec->catfile($path, $name."_tirs.gff3");
    my @tirv_opts = qw(-seqids -md5 -mintirlen -mintirdist -index -hmms);

    my @tirv_args = ("no","yes","27","100",$index,$hmmdb);
    @tirv_cmd{@tirv_opts} = @tirv_args;
    
    $self->run_tirvish(\%tirv_cmd, $gff);

    my $filtered = $self->_filter_tir_gff($gff);
    my $fixgff   = $self->_fix_tir_gff($filtered, $genome);
    #my $gff_sort = $self->sort_gff($fixgff); # no need to call sort, we can sort it

    $self->clean_index if $self->clean;
    
    return $fixgff;
}

sub _filter_tir_gff {
    my $self = shift;
    my ($gff) = @_;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_filtered.gff3");
    open my $out, '>', $outfile;

    my ($header, $features) = $self->collect_gff_features($gff);
    say $out $header;

    my @rt = qw(rve rvt rvp gag chromo);
    
    for my $tir (nsort keys %$features) {
	for my $feat (@{$features->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    if ($feats[2] eq 'protein_match') {
		my ($type, $pdom) = ($feats[8] =~ /(name) ("?\w+"?)/);
		$pdom =~ s/"//g;
		my $dom = lc($pdom);
		if ($dom =~ /rve|rvt|rvp|gag|chromo|rnase|athila|zf/i) {
		    delete $features->{$tir};
		}
	    }
	}
    }

    for my $tir (keys %$features) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s);
	my $region = @{$features->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	for my $feat (@{$features->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =~ s/\s\;\s/\;/g;
	    $feats[8] =~ s/\s+$//;
	    $feats[8] =~ s/\"//g;
	    $feats[8] =~ s/(\;\w+)\s/$1=/g;
	    $feats[8] =~ s/\s;/;/;
	    $feats[8] =~ s/^(\w+)\s/$1=/;
	    say $out join "\t", @feats;
	}
    }
    return $outfile;
}

sub _fix_tir_gff {
    my $self = shift;
    my ($gff, $genome) = @_;

    my ($idmap, $seqlen) = $self->_get_seq_len($genome);
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_id.gff3");
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";;

    my (%md5, %features);
    while (<$in>) {
        chomp;
        if (/^##sequence-region/) {
            ##sequence-region   md5:faf622ff75def8ca86a6ac12894a87dc:seq0 1 119871
            my @regs = split /\s+/;
            my ($md5, $hash, $id) = split /\:/, $regs[1];
	    if (exists $seqlen->{$regs[3]}) {
		$md5{$id} = $seqlen->{$regs[3]};
	    }
        }
        if (/^md5/) {
            my @f = split /\t/;
            my $source = shift @f;
            my ($md5, $hash, $id) = split /\:/, $source;
            $features{ $md5{$id} }{$f[3]} = join "||", @f;
        }
    }
    close $in;

    say $out "##gff-version 3";
    for my $id (nsort keys %md5) {
	say $out join q{ }, "##sequence-region", $md5{$id}, "1", $idmap->{ $md5{$id} };
    }

    for my $chr (nsort keys %features) {
        for my $start (sort { $a <=> $b } keys %{$features{$chr}}) {       
            my @feats = split /\|\|/, $features{$chr}{$start};
            say $out join "\t", $chr, @feats;
        }
    }
    close $out;

    move $outfile, $gff or die "Move failed: $!";
    return $outfile;
}

sub _get_seq_len {
    my $self = shift;
    my ($genome) = @_;
    
    my %ids;

    my $seq_in = Bio::SeqIO->new(-file => $genome, -format => 'fasta');

    while ( my $seq = $seq_in->next_seq() ) {
        my $id = $seq->id;
	$ids{$id} = $seq->length;
    }       
    
    my %len = reverse %ids;
    return (\%ids, \%len);
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

    perldoc Tephra::TIR::TIRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
