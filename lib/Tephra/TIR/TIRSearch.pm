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
use Tephra::Config::Exe;
use namespace::autoclean;

with 'Tephra::Role::Run::GT',
     'Tephra::Role::GFF',
     'Tephra::Role::Util';

has outfile => (
      is        => 'ro',
      isa       => 'Maybe[Str]',
      predicate => 'has_outfile',
      required  => 0,
);

=head1 NAME

Tephra::TIR::TIRSearch - Find TIR transposons in a reference genome

=head1 VERSION

Version 0.02.3

=cut

our $VERSION = '0.02.3';
$VERSION = eval $VERSION;

sub tir_search {
    my $self = shift;
    my ($index) = @_;
    
    my $genome  = $self->genome->absolute;
    my $hmmdb   = $self->hmmdb;
    my $outfile = $self->outfile;
    my (%suf_args, %tirv_cmd);

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($samtools) = @{$config}{qw(samtools)};

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $gff = defined $outfile ? $outfile : File::Spec->catfile($path, $name."_tirs.gff3");
    my $fas = $gff;
    $fas =~ s/\.gff.*/.fasta/;

    my @tirv_opts = qw(-seqids -mintirlen -mintirdist -index -hmms);

    my @tirv_args = ("yes","27","100",$index,$hmmdb);
    @tirv_cmd{@tirv_opts} = @tirv_args;
    
    $self->run_tirvish(\%tirv_cmd, $gff);
    
    my $filtered = $self->_filter_tir_gff($samtools, $gff, $fas);

    $self->clean_index if $self->clean;
    
    return $filtered;
}

sub _filter_tir_gff {
    my $self = shift;
    my ($samtools, $gff, $fas) = @_;
    my $genome = $self->genome;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_filtered.gff3");
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>>', $fas or die "\nERROR: Could not open file: $fas\n";
    
    my %tirs;
    my ($header, $features) = $self->collect_gff_features($gff);
    say $out $header;

    my $idx = $genome.'.fai';
    unless (-e $idx) {
	$self->_index_ref($samtools);
    }

    my @rt = qw(rve rvt rvp gag chromo);
    
    for my $tir (nsort keys %$features) {
	for my $feat (@{$features->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    if ($feats[2] eq 'terminal_inverted_repeat_element') {
		my ($id) = ($feats[8] =~ /ID\s+?=?(terminal_inverted_repeat_element\d+)/);
		$tirs{$feats[0]}{$id} = join "||", @feats[3..4];
	    }
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

    for my $tir (nsort keys %$features) {
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
	    if ($feats[2] eq 'terminal_inverted_repeat_element') {
		my ($elem) = ($feats[8] =~ /ID=(terminal_inverted_repeat_element\d+)/);
                my $tmp = $elem.".fasta";
                my $id  = $elem."_".$loc."_".$feats[3]."_".$feats[4];
                my $cmd = "$samtools faidx $genome $loc:$feats[3]-$feats[4] > $tmp";
                $self->run_cmd($cmd);
                my $seqio = Bio::SeqIO->new(-file => $tmp, -format => 'fasta');
                while (my $seqobj = $seqio->next_seq) {
                    my $seq = $seqobj->seq;
                    $seq =~ s/.{60}\K/\n/g;
                    say $faout join "\n", ">".$id, $seq;
                }
                unlink $tmp;
            }
	    say $out join "\t", @feats;
	}
    }
    close $out;

    return $outfile;
}

sub _index_ref {
    my $self = shift;
    my $fasta = $self->genome;
    my ($samtools) = @_;
    my $faidx_cmd = "$samtools faidx $fasta";
    $self->run_cmd($faidx_cmd);
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
