package Tephra::Classify::LTRSfams;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use IPC::System::Simple qw(capture);
use List::UtilsBy       qw(nsort_by);
use List::MoreUtils     qw(any none);
use Path::Class::File;
use Try::Tiny;
use Cwd;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has repeatdb => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has gff => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

#
# methods
#
sub find_gypsy_copia {
    my $self = shift;
    my ($features) = @_;
    my $is_gypsy = 0;
    my $is_copia = 0;
    my %gypsy;
    my %copia;

    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;

    my @cop_exp = qw(gag asp rve UBN2 UBN2_2 UBN2_3 RVT_2 RNase_H);
    my @gyp_exp = qw(gag asp RVT_1 RNash_H rve Chromo);
	
    for my $ltr (keys %$features) {
	my $region = @{$features->{$ltr}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$features->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    if ($feats[2] eq 'protein_match') {
		my ($doms) = ($feats[8] =~ /name \"?(\w+)\"?/);
		push @all_pdoms, $doms;
	    }
	}
	if (@all_pdoms) {
	    $pdom_org = join ",", @all_pdoms;
	    @all_pdoms = reverse @all_pdoms if $strand eq '-';
	    # gypsy -> gag,ap,int,rt,rh
	    # copia -> gag,ap,rt,rh,int
	    #
	    # model names: gag|retrotransposon_gag
	    # NAME  UBN2_2
	    # DESC  gag-polypeptide of LTR copia-type
	    # NAME  UBN2_3
	    # DESC  gag-polypeptide of LTR copia-type
	    # NAME  UBN2
	    # DESC  gag-polypeptide of LTR copia-type
	    # NAME  Retrotrans_gag
	    # DESC  Retrotransposon gag protein
	    
	    # NAME  gag-asp_proteas
	    # DESC  gag-polyprotein putative aspartyl protease
	    
	    # NAME  gag_pre-integrs
	    # DESC  GAG-pre-integrase domain
	    # NAME  rve
	    # DESC  Integrase core domain

	    # NAME  RVT_1
	    # DESC  Reverse transcriptase (RNA-dependent DNA polymerase)
	    # NAME  RVT_2
	    # DESC  Reverse transcriptase (RNA-dependent DNA polymerase)
	    # NAME  RVT_3
	    # DESC  Reverse transcriptase-like

	    # NAME  RNase_H
	    # DESC  RNase H

	    # NAME  Chromo
	    # DESC  Chromo (CHRromatin Organisation MOdifier) domain
	    #say STDERR join q{ }, "strand: $strand", $pdom_org;
	    my ($gyp_org, $cop_org);
	    my ($gyp_dom_ct, $cop_dom_ct) = (0, 0);
	    if (grep { /rvt_2|ubn/i && ! /rvt_1|chromo/i } @all_pdoms) { # &&
		#grep { ! /rvt_1|chromo/i } @all_pdoms) {
		#my @doms = any { $_ =~ @all_poms } @gyp_exp;	
		for my $d (@cop_exp) {
		    for my $p (@all_pdoms) {
			$cop_dom_ct++ if $p =~ /$d/i;
		    }
		}
	    }

	    if (grep { /rvt_1|chromo/i && ! /rvt_2|ubn/i } @all_pdoms) { # &&
		#none { /rvt_2|ubn/i } @all_pdoms) {
		for my $d (@gyp_exp) {
		    for my $p (@all_pdoms) {
			$gyp_dom_ct++ if $p =~ /$d/i;
		    }
		}
	    }

	    if ($gyp_dom_ct >= 1) {
		#say STDERR join q{ }, "Gypsy strand: $strand", $pdom_org;
		$gypsy{$ltr} = $features->{$ltr};
		delete $features->{$ltr};
	    }
	    elsif ($cop_dom_ct >= 1) {
		#say STDERR join q{ }, "Copia strand:$strand", $pdom_org;
		$copia{$ltr} = $features->{$ltr};
		delete $features->{$ltr};
	    }

	    $gyp_dom_ct = 0;
	    $cop_dom_ct = 0;
	    $is_gypsy   = 0;
	    $is_copia   = 0;

	    undef $pdom_org;
	    @all_pdoms = ();
	}
	#undef $pdom_org;
	#@all_pdoms = ();
    }    
    
    return (\%gypsy, \%copia);
}

sub find_unclassified {
    my $self = shift;
    my ($features) = @_;
    my $gff   = $self->gff;
    my $fasta = $self->genome;
    my %ltr_rregion_map;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfast = File::Spec->catfile($path, $name."_unclassified.fasta");

    open my $ofas, '>>', $outfast;

    for my $ltr (keys %$features) {
	my $region = @{$features->{$ltr}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];

	for my $feat (@{$features->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    if ($feats[2] eq 'LTR_retrotransposon') {
		my ($start, $end) = @feats[3..4];
		my ($elem) = ($feats[8] =~ /(LTR_retrotransposon\d+)/);
		my $id = $elem."_".$loc."_".$start."_".$end;
		$ltr_rregion_map{$id} = $ltr;
		my $tmp = $elem.".fasta";
		my $cmd = "samtools faidx $fasta $loc:$start-$end > $tmp";
		$self->run_cmd($cmd);
		my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
		while (my $seqobj = $seqio->next_seq) {
		    my $seq = $seqobj->seq;
		    $seq =~ s/.{60}\K/\n/g;
		    say $ofas join "\n", ">".$id, $seq;
		}   
		unlink $tmp;
	    }
	}
    }
    close $ofas;

    return ($outfast, \%ltr_rregion_map);
}

sub search_unclassified {
    my $self = shift;
    my ($unc_fas) = @_;
    my $repeatdb = $self->repeatdb;
    my $blastdb  = $self->_make_blastdb($repeatdb);

    my ($bname, $bpath, $bsuffix) = fileparse($unc_fas, qr/\.[^.]*/);
    my ($fname, $fpath, $fsuffix) = fileparse($unc_fas, qr/\.[^.]*/);
    my $outfile = $fname."_".$bname.".bln";

    my $blastcmd = "blastn -dust no -query $unc_fas -evalue 10 -db $blastdb ".
	"-outfmt 6 -num_threads 12 | sort -nrk12,12 | sort -k1,1 -u > $outfile";

    $self->run_cmd($blastcmd);

    unlink glob("$blastdb*");

    return $outfile;
}

sub annotate_unclassified {
    my $self = shift;
    my ($blast_out, $gypsy, $copia, $features, $ltr_rregion_map) = @_;
    open my $in, '<', $blast_out;
    my (%gypsy_re, %copia_re);

    while (<$in>) {
	chomp;
	my @f = split /\t/;
 
	if ($f[2] >= 80 && $f[3] >= 80) { #make length and pid thresholds options
	    my ($family) = ($f[1] =~ /(^RL[GCX][_-][a-zA-Z]*\d*?)/);
	    if ($family =~ /^RLG/) {
		$gypsy->{ $ltr_rregion_map->{$f[0]} } = $features->{ $ltr_rregion_map->{$f[0]} };
		delete $features->{ $ltr_rregion_map->{$f[0]} };
	    }
	    elsif ($family =~ /^RLC/) {
		$copia->{ $ltr_rregion_map->{$f[0]} } = $features->{ $ltr_rregion_map->{$f[0]} };
		delete $features->{ $ltr_rregion_map->{$f[0]} };
	    }
	}
    }
    close $in;
    unlink $blast_out;
}

sub write_gypsy {
    my $self = shift;
    my ($gypsy, $header) = @_;
    my $gff = $self->gff;

    my @lengths;
    my $gyp_feats;
    my %pdom_index;
    my $pdom_org;
    my $has_pdoms  = 0;
    my $pdoms      = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_gypsy.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_gypsy_domain_org.tsv");
    open my $out, '>>', $outfile;
    open my $domf, '>>', $domoutfile;
    say $out $header;

    my @all_pdoms;
    for my $ltr (nsort_by { m/repeat_region\d+\.(\d+)\.\d+/ and $1 } keys %$gypsy) {
	my ($rreg, $s, $e) = split /\./, $ltr;
	my $region = @{$gypsy->{$ltr}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];

	for my $feat (@{$gypsy->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] = $self->_format_attribute($feats[8]);
	    #if ($feats[8] =~ /Parent=repeat_region\d+(_\d+)/i) {
		#$feats[8] =~ s/$1//g;
	    #}
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=(\w+)/);
		push @all_pdoms, $doms;
	    }
	    if ($feats[2] =~ /LTR_retrotransposon/) {
		my $ltrlen = $feats[4] - $feats[3] + 1;
		push @lengths, $ltrlen;
	    }
	    $gyp_feats .= join "\t", @feats, "\n";
	}
	chomp $gyp_feats;
	say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg";
	say $out $gyp_feats;

	$pdom_org = join ",", @all_pdoms;
	$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	$pdoms++ if $has_pdoms;
	undef $gyp_feats;
	undef $pdom_org;
	@all_pdoms = ();
	$has_pdoms  = 0;
    }
    close $out;

    #dd \%pdom_index and exit;
    my %tot_dom_ct;
    say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
    for my $strand (keys %pdom_index) {
	for my $org (keys %{$pdom_index{$strand}}) {
	    $tot_dom_ct{$org} += $pdom_index{$strand}{$org};
	    say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	}
    }

    say $domf "==========";
    say $domf join "\t", "Domain_organization", "Domain_count";
    for my $domorg (keys %tot_dom_ct) {
	say $domf join "\t", $domorg, $tot_dom_ct{$domorg};
    }
    close $domf;
    
    say STDERR join "\t", "gypsy_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;

    return $outfile;
}

sub write_copia {
    my $self = shift;
    my ($copia, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $cop_feats;
    my $has_pdoms  = 0;
    my $pdoms      = 0;
    my $pdom_org;
    my @all_pdoms;
    my %pdom_index;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_copia.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_copia_domain_org.tsv");
    open my $out, '>>', $outfile;
    open my $domf, '>>', $domoutfile;
    say $out $header;

    for my $ltr (nsort_by { m/repeat_region\d+\.(\d+)\.\d+/ and $1 } keys %$copia) {
	my ($rreg, $s, $e) = split /\./, $ltr;
	my $region = @{$copia->{$ltr}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$copia->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] = $self->_format_attribute($feats[8]);
	    #if ($feats[8] =~ /Parent=repeat_region\d+(_\d+)/i) {
		#$feats[8] =~ s/$1//g;
	    #}
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=(\w+)/);
		push @all_pdoms, $doms;
		#push @{$pdom_index{$strand}}, $pdom_org;
	    }
	    if ($feats[2] =~ /LTR_retrotransposon/) {
		my $ltrlen = $feats[4] - $feats[3] + 1;
		push @lengths, $ltrlen;
	    }
	    $cop_feats .= join "\t", @feats, "\n";
	}
	chomp $cop_feats;
	say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg";
	say $out $cop_feats;

	$pdom_org = join ",", @all_pdoms;
	$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	$pdoms++ if $has_pdoms;
	undef $cop_feats;
	undef $pdom_org;
	@all_pdoms = ();
	$has_pdoms  = 0;
    }
    close $out;

    my %tot_dom_ct;
    say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
    for my $strand (keys %pdom_index) {
	for my $org (keys %{$pdom_index{$strand}}) {
	    $tot_dom_ct{$org} += $pdom_index{$strand}{$org};
	    say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	}
    }

    say $domf "==========";
    say $domf join "\t", "Domain_organization", "Domain_count";
    for my $domorg (keys %tot_dom_ct) {
	say $domf join "\t", $domorg, $tot_dom_ct{$domorg};
    }
    close $domf;
    
    say STDERR join "\t", "copia_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;

    return $outfile;
}

sub write_unclassified {
    my $self = shift;
    my ($features, $header) = @_;
    my $gff = $self->gff;

    my %pdom_index;
    my @all_pdoms;
    my $pdom_org;
    my @lengths;
    my $unc_feats;
    my $has_pdoms = 0;
    my $pdoms = 0;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_unclassified.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_unclassified_domain_org.tsv");
    open my $out, '>>', $outfile;
    open my $domf, '>>', $domoutfile;
    say $out $header;

    for my $ltr (nsort_by { m/repeat_region\d+\.(\d+)\.\d+/ and $1 } keys %$features) {
	my ($rreg, $s, $e, $l) = split /\./, $ltr;
	my $region = @{$features->{$ltr}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$features->{$ltr}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] = $self->_format_attribute($feats[8]);
	    #if ($feats[8] =~ /Parent=repeat_region\d+(_\d+)/i) {
		#$feats[8] =~ s/$1//g;
	    #}
	    if ($feats[2] =~ /protein_match/) {
		my ($doms) = ($feats[8] =~ /name=(\w+)/);
		push @all_pdoms, $doms;
		$has_pdoms = 1;
	    }
	    if ($feats[2] =~ /LTR_retrotransposon/) {
		my $ltrlen = $feats[4] - $feats[3] + 1;
		push @lengths, $ltrlen;
	    }
	    $unc_feats .= join "\t", @feats, "\n";
	}
	chomp $unc_feats;
	say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg";
	say $out $unc_feats;

	$pdom_org = join ",", @all_pdoms;
	$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	$pdoms++ if $has_pdoms;
	@all_pdoms = ();

	$has_pdoms = 0;
	undef $unc_feats;
	undef $pdom_org;
    }
    close $out;

    my %tot_dom_ct;
    say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
    for my $strand (keys %pdom_index) {
	for my $org (keys %{$pdom_index{$strand}}) {
	    $tot_dom_ct{$org} += $pdom_index{$strand}{$org};
	    say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	}
    }

    say $domf "==========";
    say $domf join "\t", "Domain_organization", "Domain_count";
    for my $domorg (keys %tot_dom_ct) {
	say $domf join "\t", $domorg, $tot_dom_ct{$domorg};
    }
    close $domf;
    
    say STDERR join "\t", "unclassified_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;
    say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;

    return $outfile;
}

sub _make_blastdb {
    my $self = shift;
    my ($db_fas) = @_;
    my ($dbname, $dbpath, $dbsuffix) = fileparse($db_fas, qr/\.[^.]*/);

    my $db = $dbname."_blastdb";
    my $dir = getcwd();
    my $db_path = Path::Class::File->new($dir, $db);
    unlink $db_path if -e $db_path;

    try {
	my @makedbout = capture([0..5],"makeblastdb -in $db_fas -dbtype nucl -title $db -out $db_path 2>&1 > /dev/null");
    }
    catch {
	say STDERR "Unable to make blast database. Here is the exception: $_.";
	say STDERR "Ensure you have removed non-literal characters (i.e., "*" or "-") in your repeat database file.";
	say STDERR "These cause problems with BLAST+. Exiting.";
	exit(1);
    };

    return $db_path;
}

sub _format_attribute {
    my $self = shift;
    my ($str) = @_;

    $str =~ s/\s\;\s/\;/g;
    $str =~ s/\s+/=/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    
    return $str;
}

__PACKAGE__->meta->make_immutable;

1;
