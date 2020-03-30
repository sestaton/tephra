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
use List::Util      qw(min max);
use List::MoreUtils qw(uniq);
use File::Path      qw(make_path);
use File::Temp      qw(tempfile);
use Cwd             qw(abs_path);         
use Carp 'croak';
use Tephra::Annotation::MakeExemplars;
use Tephra::Annotation::Util;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util',
     'Tephra::Role::Logger',
     'Tephra::Classify::Fams::Cluster',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Classify::Any - Classify any transposons into families based on similarity

=head1 VERSION

Version 0.12.5

=cut

our $VERSION = '0.12.5';
$VERSION = eval $VERSION;

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

has type => (
      is       => 'ro',
      isa      => 'Str',
      required => 1,
      default  => 'non-LTR',
);

#
# methods
#
sub process_blast_args {
    my $self = shift;
    my ($famfile) = @_;
    my $threads = $self->threads;

    my ($query, $db) = ($famfile, $famfile);
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
    elsif (+($threads-1) % 3 == 0) {
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
	#DHH_helitron1_singleton_family0_Contig57_HLAC-254L24_106214_107555
	if ($hitlen >= $blast_hlen && $hitlen >= ($minlen * $perc_cov) && $pid >= $blast_hpid) {
	    unless (exists $seen{$queryid} || exists $seen{$hitid}) {
		push @{$matches{$hitid}}, $queryid;
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
    my ($tefas, $matches, $sf_elem_map, $famct, $singct, $famtot) = @_;
    my $tetype = $self->type;

    ($famct, $singct, $famtot) = (0, 0, 0);
    #say "MATCHES ===>";
    #dd $matches;
    #say "SF_ELEM_MAP ===>";
    #dd $sf_elem_map;
    my ($name, $path, $suffix) = fileparse($tefas, qr/\.[^.]*/);
    my $dir = basename($path);

    my $seqstore = $self->store_seq($tefas);
    #my $ct0 = keys %$seqstore;
    #say "SEQSTORE ===> $ct0";
    #dd $seqstore;
    #say join "\n", keys %$seqstore;
    my $elemct = (keys %$seqstore);
	
    my $fidx = $famct;
    my $sidx = $singct;
    my $tomerge = 0;
    my (%fastas, %annot_ids, %sfmap);
    my @seen;

    my $util = Tephra::Annotation::Util->new;

    #non_LTR_retrotransposon0  => ["non_LTR_retrotransposon51_2RHet_3818_6387"],  
    if (%$matches) {
	for my $str (reverse sort { @{$matches->{$a}} <=> @{$matches->{$b}} } keys %$matches) {
	    my ($famid) = ($str =~ /(helitron\d+|non_LTR_retrotransposon\d+)/);
	    
	    unless (defined $famid) {
		say STDERR "\n[ERROR]: Could not get element ID for '$str'\n";
		$famid = $tetype;
	    }
	    
	    # sf_elem_map: non_LTR_retrotransposon99  => "RIJ",
	    my $sfcode = $sf_elem_map->{$str};
	    
	    my $famfile = $sfcode."_family$fidx".".fasta";
	    my $outfile = File::Spec->catfile( abs_path($path), $famfile );
	    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	    
	    #$famtot++;
	    $annot_ids{$str} = $sfcode."_family$fidx";
	    $self->write_element_to_family($str, $seqstore, $out, 0, $fidx, $sfcode);
	    $famtot++;
	    push @seen, $str;
	    
	    for my $elem (@{$matches->{$str}}) {
		if (exists $seqstore->{$elem}) {
		    $famtot++;
		    $annot_ids{$elem} = $sfcode."_family$fidx";
		    $self->write_element_to_family($elem, $seqstore, $out, 0, $fidx, $sfcode);
		    push @seen, $elem;
		}
		else {
		    croak "\n[ERROR]: $elem not found in store. Exiting.";
		}
	    }
	    close $out;
	    $fastas{$outfile} = 1;
	    $fidx++;
	}
    }

    if (@seen) { 
	@seen = uniq(@seen);
	#say "SEEN ===> ".scalar(@seen);
	delete $seqstore->{$_} for @seen;
    }
    #my $ct1 = keys %$seqstore;
    #say "SEQSTORE2 ===> $ct1";
    #dd $seqstore;
    #say join "\n", keys %$seqstore;

    if (%$seqstore) {
	my $famxfile = $tetype.'_singleton_families_XXXX';
	my ($outx, $xoutfile) = tempfile( TEMPLATE => $famxfile, DIR => $path, UNLINK => 0, SUFFIX => '.fasta' );

	for my $selem (nsort keys %$seqstore) {
	    my $sfcode = $sf_elem_map->{$selem};

	    $self->write_element_to_family($selem, $seqstore, $outx, 1, $sidx, $sfcode);
	    $annot_ids{$selem} = $sfcode."_singleton_family$sidx";
	    $sidx++;
	}
	close $outx;
	$fastas{$xoutfile} = 1;

    }
    #say "ANNOT_IDS ===>";    
    #dd \%annot_ids;

    return (\%fastas, \%annot_ids,
	    { total_elements    => $elemct,
	      families          => $fidx,
	      total_in_families => $famtot,
	      singletons        => $sidx });
}

sub write_element_to_family {
    my $self = shift;
    my ($elem, $seqstore, $outfh, $is_singleton, $idx, $sfcode) = @_;

    $seqstore->{$elem} =~ s/.{60}\K/\n/g;
    my ($id) = ($elem =~ /(helitron\d+|non_LTR_retrotransposon\d+)_/);
    my ($start, $stop) = ($elem =~ /(\d+)_(\d+)$/);
    my $chr = $elem;
    $chr =~ s/${id}_//;
    $chr =~ s/_$start.*//;
    my $fasid = $is_singleton ? join "_", $sfcode, "singleton_family$idx", $id, $chr, $start, $stop 
	: join "_", $sfcode, "family$idx", $id, $chr, $start, $stop;

    say $outfh join "\n", ">$fasid", $seqstore->{$elem};

    return;
}

sub combine_families {
    my ($self) = shift;
    my ($family_map, $outfile) = @_;
    
    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    
    my $ct = 0;
    for my $fam (nsort keys %$family_map) {
	for my $file (nsort keys %{$family_map->{$fam}{FAMS}}) {
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
	}
    }
    close $out;
    
    return $ct;
}

sub annotate_gff {
    my $self = shift;
    my ($family_map, $ingff, $sf_elem_map) = @_;
    my $outgff = $self->gff;

    #say "===> FAMILY_MAP";
    #dd $family_map;
    #say "===> SF_ELEM_MAP";
    #dd $sf_elem_map;

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
		#my ($fam, $id) = ($f[8] =~ /ID=((\w{3}_)?helitron\d+|(?:\w{3}_)?non_LTR_retrotransposon\d+);/);
		#my ($fam) = ($f[8] =~ /ID=(\w{3})?_helitron\d+|(\w{3})?_non_LTR_retrotransposon\d+/);
		my $id = $f[8];
		$id =~ s/ID=//;
		$id =~ s/\;.*//;
		my $key  = join "_", $id, $f[0], $f[3], $f[4];
		my $sfamily = $sf_elem_map->{$key};
		if (exists $family_map->{$sfamily}{IDS}{$key}) {
		    my $family = $family_map->{$sfamily}{IDS}{$key};
		    $f[8] =~ s/ID=$id\;/ID=$id;family=$family;/;
		    say $out join "\t", @f;
		}
		else {
		    say STDERR "\n[WARNING]: '$id' not found in Tephra::Classify::Any::annotate_gff. This is a bug. Please report it.\n";
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
