package Tephra::Classify::Any;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Sort::Naturally;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use List::Util  qw(min max);
use File::Path  qw(make_path);
use Cwd         qw(abs_path);         
#use Log::Any    qw($log);
use Carp 'croak';
use Tephra::Annotation::MakeExemplars;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util',
     'Tephra::Role::Logger',
     'Tephra::Classify::Fams::Cluster',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Classify::Any - Classify any transposons into families based on similarity

=head1 VERSION

Version 0.12.0

=cut

our $VERSION = '0.12.0';
$VERSION = eval $VERSION;

has fasta => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has outdir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 1,
    coerce   => 1,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

has gff => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

#has type => (
#      is       => 'ro',
#      isa      => 'Str',
#      required => 1,
#      default  => 'LTR',
#);

#
# methods
#
sub process_blast_args {
    my $self = shift;
    my ($obj) = @_;

    my $threads = $self->threads;
    my ($query, $db) = ($self->fasta, $self->fasta);
    unless (-s $query) {
	say STDERR "\n[WARNING]: Input FASTA is empty so the BLAST analysis will be skipped.\n";
	#unlink $query;
	return undef;
    }

    my (@fams, %exemplars);

    my $thr;
    if ($threads % 3 == 0) {
        $thr = sprintf("%.0f",$threads/2);
    }
    elsif ($threads-1 % 3 == 0) {
        $thr = sprintf("%.0f",$threads-1/3);
    }
    else {
        $thr = 1;
    }

    my $blastdb = $self->make_blastdb($db);
    my ($dbname, $dbpath, $dbsuffix) = fileparse($blastdb, qr/\.[^.]*/);
    my ($qname, $qpath, $qsuffix) = fileparse($query, qr/\.[^.]*/);
    my $report = File::Spec->catfile( abs_path($qpath), $qname."_$dbname".'.bln' );

    my $blast_report = $self->run_blast({ query   => $query, 
					  db      => $blastdb, 
					  threads => $thr, 
					  outfile => $report, 
					  sort    => 'bitscore',
					  evalue  => 1e-10 });

    my @dbfiles = glob "$blastdb*";
    unlink @dbfiles;
    #unlink $query;

    return $blast_report;
}

sub parse_blast {
    my $self = shift;
    my ($blast_report) = @_;
    return undef unless defined $blast_report;

    my $blast_hpid = 80; #$self->blast_hit_pid;
    my $blast_hcov = 50; #$self->blast_hit_cov;
    my $blast_hlen = 80; #$self->blast_hit_len;
    my $perc_cov   = sprintf("%.2f", $blast_hcov/100);

    my (%matches, %seen);
    open my $in, '<', $blast_report or die "\n[ERROR]: Could not open file: $blast_report\n";
    while (my $line = <$in>) {
	chomp $line;
	my ($queryid, $hitid, $pid, $hitlen, $mmatchct, $gapct, 
	    $qhit_start, $qhit_end, $hhit_start, $hhit_end, $evalue, $score) = split /\t/, $line;
	next if $queryid eq $hitid;
	my ($qstart, $qend) = ($queryid =~ /(\d+)-?_?(\d+)$/);
	my ($hstart, $hend) = ($hitid =~ /(\d+)-?_?(\d+)$/);
	my $qlen = $qend - $qstart + 1;
	my $hlen = $hend - $hstart + 1;
	my $minlen = min($qlen, $hlen); # we want to measure the coverage of the smaller element
	#my ($coords) = ($queryid =~ /_(\d+_\d+)$/);
        #$queryid =~ s/_$coords//;
	#DHH_helitron1_singleton_family0_Contig57_HLAC-254L24_106214_107555
	my ($elem) = ($hitid =~ /(non_LTR_retrotransposon\d+|helitron\d+)_/);
	#my ($family, $elem) = ($hitid =~ /^(\w{3})_(non_LTR_retrotransposon\d+|helitron\d+)_/);
	if ($hitlen >= $blast_hlen && $hitlen >= ($minlen * $perc_cov) && $pid >= $blast_hpid) {
	    unless (exists $seen{$queryid}) {
		push @{$matches{$elem}}, $queryid;
		$seen{$queryid} = 1;
	    }
	}
    }
    close $in;
    unlink $blast_report;

    return \%matches;
}

sub write_families {
    my $self = shift;
    my ($tefas, $matches, $sf_elem_map, $sf) = @_;

    my ($name, $path, $suffix) = fileparse($tefas, qr/\.[^.]*/);
    my $dir  = basename($path);
    #my ($sf) = ($dir =~ /_(\w+)$/);
    #my ($sf) = ($dir =~ /_((?:\w+\d+-)?\w+)$/);
    #unless (defined $sf) {
        #say STDERR "\n[ERROR]: Can't get sf from $dir $.";
    #}

    my $seqstore = $self->_store_seq($tefas);
    #dd $seqstore;
    my $elemct = (keys %$seqstore);
	
    my ($idx, $famtot, $tomerge) = (0, 0, 0);
    my (%fastas, %annot_ids, %sfmap);

    for my $str (reverse sort { @{$matches->{$a}} <=> @{$matches->{$b}} } keys %$matches) {
	my $famfile = $sf."_family$idx".".fasta";
	my $outfile = File::Spec->catfile( abs_path($path), $famfile );
	open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	for my $elem (@{$matches->{$str}}) {
	    if (exists $seqstore->{$elem}) {
		$famtot++;
		$seqstore->{$elem} =~ s/.{60}\K/\n/g;
		$annot_ids{$elem} = "family$idx";
		my ($id) = ($elem =~ /(helitron\d+|non_LTR_retrotransposon\d+)_/);
		my ($start, $stop) = ($elem =~ /(\d+)_(\d+)$/);
		my $chr = $elem;
		$chr =~ s/${id}_//;
		$chr =~ s/_$start.*//;
		my $sfcode = $sf_elem_map->{$id};
		#my $sfcode = $sf_elem_map->{$elem};
		#say $out join "\n", ">$elem"."_family$idx", $seqstore->{$elem};
		say $out join "\n", ">$sfcode"."_family$idx"."_$id"."_$chr"."_$start"."_$stop", $seqstore->{$elem};
		delete $seqstore->{$elem};
	    }
	    else {
		croak "\n[ERROR]: $elem not found in store. Exiting.";
	    }
	}
	close $out;
	$idx++;
	$fastas{$outfile} = 1;
    }
    my $famct = $idx;
    $idx = 0;

    #dd $sf_elem_map;
    if (%$seqstore) {
	my $famxfile = $sf.'_singleton_families.fasta';
	my $xoutfile = File::Spec->catfile( abs_path($path), $famxfile );
	open my $outx, '>', $xoutfile or die "\n[ERROR]: Could not open file: $xoutfile\n";
	for my $k (nsort keys %$seqstore) {
	    $seqstore->{$k} =~ s/.{60}\K/\n/g;
	    my $chr = $k;
	    my ($id) = ($k =~ /(helitron\d+|non_LTR_retrotransposon\d+)_/); 
	    #my ($sfcode, $id) = ($k =~ /(\w{3})_(helitron\d+|non_LTR_retrotransposon\d+)_/);
	    my ($start, $stop) = ($k =~ /(\d+)_(\d+)$/);
	    $chr =~ s/${id}_//;
	    $chr =~ s/_$start.*//;
	    my $sfcode = $sf_elem_map->{$id};
	    #say $outx join "\n", ">$id"."_singleton_family$idx"."_$chr"."_$start"."_$stop", $seqstore->{$k};
	    say $outx join "\n", ">$sfcode"."_singleton_family$idx"."_$id"."_$chr"."_$start"."_$stop", $seqstore->{$k};
	    $annot_ids{$k} = "singleton_family$idx";
	    $idx++;
	}
	close $outx;
	$fastas{$xoutfile} = 1;
    }
    my $singct = $idx;
    
    return (\%fastas, \%annot_ids, \%sfmap,
	    { total_elements    => $elemct,
	      families          => $famct,
	      total_in_families => $famtot,
	      singletons        => $singct });
}

sub combine_families {
    my ($self) = shift;
    my ($outfiles, $outfile) = @_; #$sf_elem_map) = @_;
    
    #dd $sf_elem_map;
    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";

    my $ct = 0;
    for my $file (nsort keys %$outfiles) {
	unlink $file && next unless -s $file;
	my $kseq = Bio::DB::HTS::Kseq->new($file);
	my $iter = $kseq->iterator();
	while (my $seqobj = $iter->next_seq) {
	    my $id  = $seqobj->name;
	    my $seq = $seqobj->seq;
	    $seq =~ s/.{60}\K/\n/g;
	    #say $out join "\n", ">$sf_elem_map->{$key}"."_$id", $seq;
	    say $out join "\n", ">$id", $seq;
	    $ct++;
	}
	#unlink $file;
    }
    close $out;

    return $ct;
}

sub annotate_gff {
    my $self = shift;
    my ($annot_ids, $ingff, $sf_elem_map) = @_;
    my $outgff = $self->gff;

    open my $in, '<', $ingff or die "\n[ERROR]: Could not open file: $ingff\n";
    open my $out, '>', $outgff or die "\n[ERROR]: Could not open file: $outgff\n";

    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^#/) {
	    say $out $line;
	}
	else {
	    my @f = split /\t/, $line;
	    if ($f[2] =~ /helitron|non_LTR_retrotransposon/) {
		my ($id) = ($f[8] =~ /ID=((?:\w{3}_)?helitron\d+|(?:\w{3}_)?non_LTR_retrotransposon\d+);/);
		my $key  = join "_", $id, $f[0], $f[3], $f[4];
		if (exists $annot_ids->{$key}) {
		    my $family = $annot_ids->{$key};
		    my $sfamily = $sf_elem_map->{$id};
		    #my $fid = $sfamily."_$id";
		    my $fid = $sfamily."_$family";
		    $f[8] =~ s/ID=$id\;/ID=$id;family=$fid;/;
		    say $out join "\t", @f;
		}
		else {
		    say $out join "\t", @f;
		}
	    }
	    else {
		say $out join "\t", @f;
	    }
	}
    }
    close $in;
    close $out;

    return;
}

sub _store_seq {
    my $self = shift;
    my ($file) = @_;

    my %hash;
    my $kseq = Bio::DB::HTS::Kseq->new($file);
    my $iter = $kseq->iterator();
    while (my $seqobj = $iter->next_seq) {
	my $id = $seqobj->name;
	my $seq = $seqobj->seq;
	$hash{$id} = $seq;
    }

    return \%hash;
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

    perldoc Tephra::Classify::LTRFams


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 
	
__PACKAGE__->meta->make_immutable;

1;
