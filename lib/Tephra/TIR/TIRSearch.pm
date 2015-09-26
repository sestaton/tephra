package Tephra::TIR::TIRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use IPC::System::Simple qw(system EXIT_ANY);
use Sort::Naturally;
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;
use namespace::autoclean;

with 'Tephra::Role::Run::GT',
     'Tephra::Role::GFF';

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
    my $gff_sort = $self->sort_gff($fixgff);

    $self->clean_index if $self->clean;
    
    return $gff_sort;
}

sub _filter_tir_gff {
    my $self = shift;
    my ($gff) = @_;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_filtered.gff3");
    open my $out, '>', $outfile;

    my ($header, $features) = $self->collect_gff_features($gff);
    say $out $header;

    my $ltrrt = 0;
    my @rt = qw(rve rvt rvp gag chromo);
    
    for my $tir (nsort keys %$features) {
	#$tirct++;
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
	#push @lengths, $len;
	my $region = @{$features->{$tir}}[0];
	my ($loc, $source) = (split /\|\|/, $region)[0..1];
	say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	for my $feat (@{$features->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    #$feats[8] =~ s/\s\;\s/\;/g;
	    #$feats[8] =~ s/\s+/=/g;
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

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_id.gff3");
    open my $out, '>', $outfile;
    open my $in, '<', $gff;

    my $seqio = Bio::SeqIO->new( -file => $genome, -format => 'fasta' );

    my (%md5, %name_map, $seqid, $seqlen, $first_id);
    my $write_id = 1;
    while (<$in>) {
	chomp;
	if (/^##gff/) {
	    say $out $_;
	}
	elsif (/^##sequence-region/) {
	    ##sequence-region   md5:faf622ff75def8ca86a6ac12894a87dc:seq0 1 119871
	    my @regs = split /\s+/;
	    $md5{$regs[1]} = $regs[3];
	}
	elsif (/^#[^#]/) {
	    say STDERR $_;
	    s/#//g;
	    my $len;
	    my @name = split /\s+/;
	    if ($name[1] =~ /len=(\d+)/) {
		$len = $1;
	    }
	    my ($key) = grep { $md5{$_} eq $len } keys %md5;
	    $name_map{$key} = $name[0];
	    say $out join q{ }, "##sequence-region", $name[0], "1", $len;
	}
	elsif (/^md5/ && (keys %md5) == 1 && $write_id) {
	    my @f = split /\t/;
	    my $len; for my $k (keys %md5) { $len = $md5{$k} }
	    while (my $seqobj = $seqio->next_seq) {
		my $seq = $seqobj->seq;
		my $seql = length($seq);
		if ($len == $seql) {
		    $seqlen = $seql;
		    $seqid  = $seqobj->id;
		}
	    }
	    say $out join q{ }, "##sequence-region", $seqid, "1", $seqlen;
	    say $out join "\t", $seqid, @f[1..$#f];
	    $write_id = 0;
	    
	}
	elsif (/^md5/ && (keys %md5) == 1 && !$write_id) {
	    my @f = split /\t/;
	    say $out join "\t", $seqid, @f[1..$#f];
	}
    }
    close $in;
    close $out;

    return $outfile;
}

__PACKAGE__->meta->make_immutable;

1;
