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
use Try::Tiny;
use Tephra::Annotation::MakeExemplars;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util',
     'Tephra::Classify::Fams::Cluster',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Classify::Any - Classify any transposons into families based on similarity

=head1 VERSION

Version 0.08.1

=cut

our $VERSION = '0.08.1';
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
    unless (-s $query && -s $db) {
	unlink $query unless -s $query;
	unlink $db unless -s $db;
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
    #unlink $query, $db;

    return $blast_report;
}

sub parse_blast {
    my $self = shift;
    my ($blast_report) = @_;
    return undef unless defined $blast_report;

    my $blast_hpid = $self->blast_hit_pid;
    my $blast_hcov = $self->blast_hit_cov;
    my $blast_hlen = $self->blast_hit_len;
    my $perc_cov   = sprintf("%.2f", $blast_hcov/100);

    my (%matches, %seen);
    open my $in, '<', $blast_report or die "\nERROR: Could not open file: $blast_report\n";
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
	my ($coords) = ($queryid =~ /_(\d+_\d+)$/);
        $queryid =~ s/_$coords//;
	my ($family) = ($hitid =~ /(\w{3}_(?:non_LTR_retrotransposon|helitron)\d+)_/);
	if ($hitlen >= $blast_hlen && $hitlen >= ($minlen * $perc_cov) && $pid >= $blast_hpid) {
	    unless (exists $seen{$queryid}) {
		push @{$matches{$family}}, $queryid;
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
    my ($matches, $tetype, $tefas) = @_;

    my ($name, $path, $suffix) = fileparse($tefas, qr/\.[^.]*/);
    #my $dir  = basename($cpath);
    #my ($sf) = ($dir =~ /_(\w+)$/);

    #my $sfname;
    #if ($tetype eq 'non-LTR') {
	#$sfname = 'RLG' if $sf =~ /gypsy/i;
	#$sfname = 'RLC' if $sf =~ /copia/i;
	#$sfname = 'RLX' if $sf =~ /unclassified/i;
    #}
    #if ($tetype eq 'Helitron') {
	#$sfname = 'DTA' if $sf =~ /hat/i;
	#$sfname = 'DTC' if $sf =~ /cacta/i;
	#$sfname = 'DTM' if $sf =~ /mutator/i;
	#$sfname = 'DTT' if $sf =~ /mariner/i;
	#$sfname = 'DTX' if $sf =~ /unclassified/i;
    #}

    #my @compfiles;
    #find( sub { push @compfiles, $File::Find::name if /complete.fasta$/ }, $cpath );
    #my $ltrfas = shift @compfiles;
    my $seqstore = $self->_store_seq($tefas);
    my $elemct = (keys %$seqstore);
	
    my ($idx, $famtot, $tomerge) = (0, 0, 0);
    my (%fastas, %annot_ids, %sfmap);

    #my $fam_id_map;
    #if (defined $matches) {
	#my $tomerge = $self->_compare_merged_nonmerged($matches, $seqstore, $unmerged_stats);
	#$fam_id_map = $tomerge == 1 ? $matches : $dom_fam_map;
    #}

    #for my $str (reverse sort { @{$fam_id_map->{$a}} <=> @{$fam_id_map->{$b}} } keys %$fam_id_map) {
    for my $str (reverse sort { @{$matches->{$a}} <=> @{$matches->{$b}} } keys %$matches) {
	my $famfile = $sf."_family$idx".".fasta";
	my $outfile = File::Spec->catfile( abs_path($path), $famfile );
	open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
	for my $elem (@{$matches->{$str}}) {
	    #my $query = $elem =~ s/\w{3}_singleton_family\d+_//r;
	    my ($sfname) = ($elem =~ /^(\w{3})_\w+/);
	    if (exists $seqstore->{$query}) {
		$famtot++;
		my $coordsh = $seqstore->{$query};
		my $coords  = (keys %$coordsh)[0];
		$seqstore->{$query}{$coords} =~ s/.{60}\K/\n/g;
		say $out join "\n", ">$sfname"."_family$idx"."_$query"."_$coords", $seqstore->{$query}{$coords};
		delete $seqstore->{$query};
		$annot_ids{$query} = $sfname."_family$idx";
		$sfmap{$sfname}++;
	    }
	    else {
		croak "\nERROR: $query not found in store. Exiting.";
	    }
	}
	close $out;
	$idx++;
	$fastas{$outfile} = 1;
    }
    my $famct = $idx;
    $idx = 0;
	
    if (%$seqstore) {
	my $famxfile = $sf.'_singleton_families.fasta';
	my $xoutfile = File::Spec->catfile( abs_path($path), $famxfile );
	open my $outx, '>', $xoutfile or die "\nERROR: Could not open file: $xoutfile\n";
	for my $k (nsort keys %$seqstore) {
	    my ($sfname) = ($k =~ /^(\w{3})_\w+/);
	    $sfmap{$sfname}++;
	    my $coordsh = $seqstore->{$k};
	    my $coords  = (keys %$coordsh)[0];
	    $seqstore->{$k}{$coords} =~ s/.{60}\K/\n/g;
	    say $outx join "\n", ">$sfname"."_singleton_family$idx"."_$k"."_$coords", $seqstore->{$k}{$coords};
	    $annot_ids{$k} = $sfname."_singleton_family$idx";
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
    my ($outfiles) = @_;
    #my $genome = $self->genome->absolute->resolve;
    #my $outdir = $self->outdir->absolute->resolve;
    my $outfile = $self->fasta;
    
    #my ($name, $path, $suffix) = fileparse($outgff, qr/\.[^.]*/);
    #my $outfile = File::Spec->catfile($path, $name.'.fasta');

    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

    my $ct = 0;
    for my $file (nsort keys %$outfiles) {
	unlink $file && next unless -s $file;
	my $kseq = Bio::DB::HTS::Kseq->new($file);
	my $iter = $kseq->iterator();
	while (my $seqobj = $iter->next_seq) {
	    my $id  = $seqobj->name;
	    my $seq = $seqobj->seq;
	    $seq =~ s/.{60}\K/\n/g;
	    say $out join "\n", ">$id", $seq;
	    $ct++;
	}
	#unlink $file;
    }
    close $outfile;

    return $ct;
}

sub annotate_gff {
    my $self = shift;
    my ($annot_ids, $ingff) = @_;
    #my $outdir = $self->outdir->absolute->resolve;
    my $outgff = $self->gff;

    open my $in, '<', $ingff or die "\nERROR: Could not open file: $ingff\n";
    open my $out, '>', $outgff or die "\nERROR: Could not open file: $outgff\n";

    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^#/) {
	    say $out $line;
	}
	else {
	    my @f = split /\t/, $line;
	    if ($f[2] eq 'helitron|non_LTR_retrotransposon') {
		my ($id) = ($f[8] =~ /ID=(helitron\d+|non_LTR_retrotransposon\d+);/);
		my $key  = $id."_$f[0]";
		if (exists $annot_ids->{$key}) {
		    my $family = $annot_ids->{$key};
		    $f[8] =~ s/ID=$id\;/ID=$id;family=$family;/;
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

sub _compare_merged_nonmerged {
    my $self = shift;
    my ($matches, $seqstore, $unmerged_stats) = @_;

    my $famtot = 0;

    for my $str (keys %$matches) {
	for my $elem (@{$matches->{$str}}) {
	    my $query = $elem =~ s/\w{3}_singleton_family\d+_//r;
	    if (exists $seqstore->{$query}) {
		$famtot++;
	    }
	    else {
		croak "\nERROR: $query not found in store. Exiting.";
	    }
	    
	}
    }

    my $tomerge = $unmerged_stats->{total_in_families} < $famtot ? 1 : 0;

    return $tomerge;
}

sub _store_seq {
    my $self = shift;
    my ($file) = @_;

    my %hash;
    my $kseq = Bio::DB::HTS::Kseq->new($file);
    my $iter = $kseq->iterator();
    while (my $seqobj = $iter->next_seq) {
	my $id = $seqobj->name;
	my ($coords) = ($id =~ /_(\d+_\d+)$/);
	$id =~ s/_$coords//;
	my $seq = $seqobj->seq;
	$hash{$id} = { $coords => $seq };
    }

    return \%hash;
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

    perldoc Tephra::Classify::LTRFams


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 
	
__PACKAGE__->meta->make_immutable;

1;
