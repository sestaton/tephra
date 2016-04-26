package Tephra::LTR::LTRRefine;

use 5.010;
use Moose;
use autodie;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use File::Copy          qw(move);
use IPC::System::Simple qw(system EXIT_ANY);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use Bio::SeqIO;
use Bio::Tools::GFF;
use Set::IntervalTree;
use Path::Class::File;
use Try::Tiny;
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::LTR::LTRRefine - Refine LTR predictions based on multiple criteria

=head1 VERSION

Version 0.02.6

=cut

our $VERSION = '0.02.6';
$VERSION = eval $VERSION;

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has outfile => (
      is        => 'ro',
      isa       => 'Str',
      predicate => 'has_outfile',
);

has n_threshold => (
      is       => 'ro',
      isa      => 'Num',
      required => 0,
      default  => 0.30,
);

has remove_dup_domains => (
    is        => 'ro',
    isa       => 'Bool',
    predicate => 'has_remove_dup_domains',
    lazy      => 1,
    default   => 0,
);

has remove_tnp_domains => (
    is        => 'ro',
    isa       => 'Bool',
    predicate => 'has_remove_tnp_domains',
    lazy      => 1,
    default   => 0,
);
#
# methods
#
sub collect_features {
    my $self = shift;
    my ($gff_thresh_ref) = @_;
    my ($gff, $pid_thresh) = @{$gff_thresh_ref}{qw(gff pid_threshold)};

    my %intervals;
    my $fasta = $self->genome;
    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

    my ($source, $start, $end, $length, $region, %features);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'repeat_region') {
	    my @string = split /\t/, $feature->gff_string;
	    $source = $string[0];
	    ($region) = ($string[8] =~ /ID=?\s+?(repeat_region\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $length = $end - $start + 1;
	}
	next $feature unless defined $start && defined $end;
	if ($feature->primary_tag ne 'repeat_region') {
	    if ($feature->start >= $start && $feature->end <= $end) {
		$intervals{$region."_".$pid_thresh} = join ".", $start, $end, $length;
		my $region_key = join ".", $region."_".$pid_thresh, $start, $end, $length;
		my @feats = split /\t/, $feature->gff_string;
		if ($feats[8] =~ /(repeat_region\d+)/) {
		    my $old_parent = $1;
		    my $new_parent = $old_parent."_".$pid_thresh;
		    $feats[8] =~ s/$old_parent/$new_parent/g;
		}
		push @{$features{$source}{$region_key}}, join "||", @feats;
	    }
	}
    }

    my ($filtered, $stats) = $self->_filter_compound_elements(\%features, $fasta);
    return ({ collected_features => $filtered, stats => $stats, intervals => \%intervals });
}

sub get_overlaps {
    my $self = shift;
    my ($feature_ref) = @_;
    my ($relaxed_features, $strict_features) = @{$feature_ref}{qw(relaxed_features strict_features)};
    my $allfeatures  = $relaxed_features->{collected_features};
    my $partfeatures = $strict_features->{collected_features}; 
    my $intervals    = $relaxed_features->{intervals};
    
    my (@best_elements, %chr_intervals);
  
    for my $source (keys %$allfeatures) {
	my $tree = Set::IntervalTree->new;
	
	for my $rregion (keys %{$allfeatures->{$source}}) {
	    my ($reg, $start, $end, $length) = split /\./, $rregion;
	    $tree->insert($reg, $start, $end);
	    $chr_intervals{$source} = $tree;
	}
    }

    for my $source (keys %$partfeatures) {
	for my $rregion (keys %{$partfeatures->{$source}}) {
	    my (%scores, %sims);

	    if (exists $chr_intervals{$source}) {
		my ($reg, $start, $end, $length) = split /\./, $rregion;
		my $res = $chr_intervals{$source}->fetch($start, $end);

		next unless defined $partfeatures->{$source}{$rregion};
		my ($score99, $sim99) = $self->_summarize_features($partfeatures->{$source}{$rregion});
		
		$scores{$source}{$rregion} = $score99;
		$sims{$source}{$rregion} = $sim99;
		
		## collect all features, then compare scores/sim...
		for my $over (@$res) {
		    my ($s, $e, $l) = split /\./, $intervals->{$over};
		    my $region_key = join ".", $over, $s, $e, $l;
		    next unless defined $allfeatures->{$source}{$region_key};
		    my ($score85, $sim85) = $self->_summarize_features($allfeatures->{$source}{$region_key});

		    $scores{$source}{$region_key} = $score85;
		    $sims{$source}{$region_key} = $sim85;
		}

		if (@$res > 0) {
		    my $best_element 
			= $self->_get_ltr_score_dups(\%scores, \%sims, $allfeatures, $partfeatures);
		    push @best_elements, $best_element;
		}
	    }
	}
    }
    
    return \@best_elements;
}

sub reduce_features {
    my $self = shift;
    my $fasta = $self->genome;
    my ($feature_ref) = @_;
    my ($relaxed_features, $strict_features, $best_elements)
	= @{$feature_ref}{qw(relaxed_features strict_features best_elements)};

    my $all_feats  = $relaxed_features->{collected_features};
    my $part_feats = $strict_features->{collected_features};
    my $all_stats  = $relaxed_features->{stats};
    my $part_stats = $strict_features->{stats};
    
    my ($all, $best, $part, $comb) = (0, 0, 0, 0);
    
    my (%best_features, %all_features, %best_stats);

    for my $stat (keys %$all_stats) {
	if (exists $part_stats->{$stat}) {
	    $best_stats{$stat} = $all_stats->{$stat} + $part_stats->{$stat};
	}
    }
    
    for my $str (@$best_elements) {
	for my $chr (keys %$str) {
	    for my $element (keys %{$str->{$chr}}) {
		$best++;
		$best_features{$chr}{$element} = $str->{$chr}{$element};
	    }
	}
    }

    for my $chromosome (keys %best_features) {
	for my $element (keys %{$best_features{$chromosome}}) {
	    if (exists $all_feats->{$chromosome}{$element}) {
		delete $all_feats->{$chromosome}{$element};
	    }
	    elsif (exists $part_feats->{$chromosome}{$element}) {
		delete $part_feats->{$chromosome}{$element};
	    }
	}
    }

    for my $source (keys %$all_feats) {
	for my $element (keys %{$all_feats->{$source}}) {
	    $all++;
	    $best_features{$source}{$element} = $all_feats->{$source}{$element};
	}
    }

    for my $source (keys %$part_feats) {
	for my $element (keys %{$part_feats->{$source}}) {
	    $part++;
	    $best_features{$source}{$element} = $part_feats->{$source}{$element};
	}
    }

    my $n_perc_filtered = 0;
    my $n_thresh = $self->n_threshold;

    for my $source (keys %best_features) {
	for my $element (keys %{$best_features{$source}}) {
	    $comb++;
	    my ($region, $start, $end, $length) = split /\./, $element;
	    my $key = join "||", $region, $start, $end;

	    my $n_perc = $self->_filterNpercent($source, $key, $fasta);
	    if ($n_perc >= $n_thresh) {
		#say STDERR "=====> Over thresh: $n_perc";
		delete $best_features{$source}{$element};
		$n_perc_filtered++;
		$comb--;
	    }
	}
    }

    $best_stats{n_perc_filtered} = $n_perc_filtered;

    print STDERR join q{ }, "Number of elements filtered by type:";
    for my $s (keys %best_stats) {
	print STDERR " $s=$best_stats{$s}";
    }

    say STDERR join q{ }, "\nNumber of elements found under what constraints:", 
        "Relaxed=$all", "Strict=$part", "Best=$best", "Combined=$comb";

    return \%best_features;
}

sub sort_features {
    my $self = shift;
    my ($feature_ref) = @_;
    my ($gff, $combined_features) = @{$feature_ref}{qw(gff combined_features)};
    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($samtools) = @{$config}{qw(samtools)};

    my $fasta = $self->genome;
    $self->_index_ref($samtools, $fasta);
    my ($outfile, $outfasta);
 
    if ($self->has_outfile) {
	$outfile = $self->outfile;
	my ($name, $path, $suffix) = fileparse($outfile, qr/\.[^.]*/);
	$outfasta = File::Spec->catfile($path, $name.".fasta");
    }
    else {
	my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
	$outfasta = File::Spec->catfile($path, $name."_combined_filtered.fasta");
	$outfile  = File::Spec->catfile($path, $name."_combined_filtered.gff3");
    }

    open my $ofas, '>>', $outfasta or die "\nERROR: Could not open file: $outfasta\n";

    my ($elem_tot, $index) = (0, 1);
    if (defined $combined_features) {
	open my $ogff, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

	my ($header, %features);
	open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";;
	while (<$in>) {
	    chomp;
	    if (/^#/) {
		$header .= $_."\n";
	    }
	    else {
		last;
	    }
	}
	close $in;
	chomp $header;
	say $ogff $header;
	
	for my $chromosome (nsort keys %$combined_features) {
	    for my $ltr (nsort_by { m/repeat_region\d+\_\d+\.(\d+)\.\d+/ and $1 }
			 keys %{$combined_features->{$chromosome}}) {
		my ($rreg, $rreg_start, $rreg_end, $rreg_length) = split /\./, $ltr;
		my $new_rreg = $rreg;
		$new_rreg =~ s/\d+.*/$index/;
		my ($first) = @{$combined_features->{$chromosome}{$ltr}}[0];
		my ($source, $strand) = (split /\|\|/, $first)[1,6];
		say $ogff join "\t", $chromosome, $source, 'repeat_region', 
	            $rreg_start, $rreg_end, '.', $strand, '.', "ID=$new_rreg";
		for my $entry (@{$combined_features->{$chromosome}{$ltr}}) {
		    my @feats = split /\|\|/, $entry;
		    $feats[8] =~ s/\s\;\s/\;/g;
		    $feats[8] =~ s/\s+$//;
		    $feats[8] =~ s/\"//g;
		    $feats[8] =~ s/(\;\w+)\s/$1=/g;
		    $feats[8] =~ s/\s;/;/;
		    $feats[8] =~ s/^(\w+)\s/$1=/;
		    $feats[8] =~ s/repeat_region\d+_\d+/$new_rreg/g;
		    $feats[8] =~ s/LTR_retrotransposon\d+/LTR_retrotransposon$index/g;
		    say $ogff join "\t", @feats;
		    
		    if ($feats[2] eq 'LTR_retrotransposon') {
			#$feats[8] .= ";SO:0000186";
			$elem_tot++;
			my ($start, $end) = @feats[3..4];
			my ($elem) = ($feats[8] =~ /(LTR_retrotransposon\d+)/);
			$elem =~ s/\d+.*//;
			$elem .= $index;
			my $id = $elem."_".$chromosome."_".$start."_".$end;
			$self->_get_ltr_range($samtools, $fasta, $elem, $id, $chromosome, $start, $end, $ofas);
		    }
		}
		$index++;
	    }
	}
	close $ogff;
	
	say STDERR "\nTotal elements written: $elem_tot";
    }
    else {
	open my $in, '<', $gff, or die "\nERROR: Could not open file: $gff\n";
	while (my $line = <$in>) {
	    chomp $line;
	    next if $line =~ /^#/;
	    my @feats = split /\t/, $line;
	    if ($feats[2] eq 'LTR_retrotransposon') {
		#$feats[8] .= ";SO:0000186";                                                            
		$elem_tot++;
		my ($chromosome, $start, $end) = @feats[0,3,4];
		my ($elem) = ($feats[8] =~ /(LTR_retrotransposon\d+)/);
		$elem =~ s/\d+.*//;
		$elem .= $index;
		my $id = $elem."_".$chromosome."_".$start."_".$end;
		$self->_get_ltr_range($samtools, $fasta, $elem, $id, $chromosome, $start, $end, $ofas);
	    }
	}
	close $in;

	move $gff, $outfile or die "Move failed: $!";
	say STDERR "\nTotal elements written: $elem_tot";
    }
    close $ofas;

    #say STDERR "\nTotal elements written: $elem_tot";
}

sub _get_ltr_range {
    my $self = shift;
    my ($samtools, $fasta, $elem, $id, $chromosome, $start, $end, $ofh) = @_;
    my $tmp = $elem.".fasta";
    my $cmd = "$samtools faidx $fasta $chromosome:$start-$end > $tmp";
    $self->run_cmd($cmd);

    my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
    while (my $seqobj = $seqio->next_seq) {
	my $seq = $seqobj->seq;
	$seq =~ s/.{60}\K/\n/g;
	say $ofh join "\n", ">".$id, $seq;
    }
    unlink $tmp;
}

sub _get_ltr_score_dups {
    my $self = shift;
    my ($scores, $sims, $allfeatures, $partfeatures) = @_;
    my %sccounts;
    my %sicounts;
    my %best_element;
    my ($score_best, $sims_best) = (0, 0);
    
    my $max_score = max(values %$scores);
    my $max_sims  = max(values %$sims);

    for my $source (sort keys %$scores) {
	for my $score_key (sort keys %{$scores->{$source}}) {
	    my $score_value = $scores->{$source}{$score_key};
	    push @{$sccounts{$score_value}}, $score_key;
	}
    }

    for my $source (sort keys %$sims) {
	for my $sim_key (sort keys %{$sims->{$source}}) {
	    my $sim_value = $sims->{$source}{$sim_key};
	    push @{$sicounts{$sim_value}}, $sim_key;
	}
    }

    my ($best_score_key, $best_sim_key);
    for my $src (keys %$scores) {
	$best_score_key = (reverse sort { $scores->{$src}{$a} <=> $scores->{$src}{$b} } keys %{$scores->{$src}})[0];
    }

    for my $src (keys %$sims) {
	$best_sim_key = (reverse sort { $sims->{$src}{$a} <=> $sims->{$src}{$b} } keys %{$sims->{$src}})[0];
    }

    for my $source (keys %$scores) {
	if (@{$sccounts{ $scores->{$source}{$best_score_key} }} == 1 &&
	    @{$sicounts{ $sims->{$source}{$best_sim_key} }} == 1  &&
	    $best_score_key eq $best_sim_key) {
	    $score_best = 1;
	    my $bscore  = $scores->{$source}{$best_score_key};
	    my $bsim    = $sims->{$source}{$best_score_key};
	    if (exists $partfeatures->{$source}{$best_score_key}) {
		$best_element{$source}{$best_score_key} = $partfeatures->{$source}{$best_score_key};
	    }
	    elsif (exists $allfeatures->{$source}{$best_score_key}) {
		$best_element{$source}{$best_score_key} = $allfeatures->{$source}{$best_score_key};
	    }
	    else {
		say "\nERROR: Something went wrong....'$best_score_key' not found in hash. This is a bug, please report it.";
		exit(1);
	    }
	}
	elsif (@{$sccounts{ $scores->{$source}{$best_score_key} }} >= 1 && 
	       @{$sicounts{ $sims->{$source}{$best_sim_key} }} >=  1) {
	    $sims_best = 1;
	    my $best   = @{$sicounts{ $sims->{$source}{$best_sim_key} }}[0];
	    my $bscore = $scores->{$source}{$best_score_key};
	    my $bsim   = $sims->{$source}{$best_score_key};
	    if (exists $partfeatures->{$source}{$best}) {
		$best_element{$source}{$best_score_key} = $partfeatures->{$source}{$best};
	    }
	    elsif (exists $allfeatures->{$source}{$best}) {
		$best_element{$source}{$best_score_key} = $allfeatures->{$source}{$best};
	    }
	    else {
		say "\nERROR: Something went wrong....'$best' not found in hash. This is a bug, please report it.";
		exit(1);
	    }
	}
    }

    if ($score_best) {
	for my $src (keys %$scores) {
	    for my $element (keys %{$scores->{$src}}) {
		if (exists $allfeatures->{$src}{$element}) {
		    delete $allfeatures->{$src}{$element};
		}
		elsif (exists $partfeatures->{$src}{$element}) {
		    delete $partfeatures->{$src}{$element};
		}
	    }
	}
    }
    elsif ($sims_best) {
	for my $src (keys %$sims) {
	    for my $element (keys %{$sims->{$src}}) {
		if (exists $allfeatures->{$src}{$element}) {
		    delete $allfeatures->{$src}{$element};
		}
		elsif (exists $partfeatures->{$src}{$element}) {
		    delete $partfeatures->{$src}{$element};
		}
	    }
	}
    }
    $score_best = 0;
    $sims_best  = 0;
    
    return \%best_element;
}

sub _summarize_features {
    my $self = shift;
    my ($feature) = @_;
    my ($three_pr_tsd, $five_pr_tsd, $ltr_sim);
    my ($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq, $tsd_ct) = (0, 0, 0, 0, 0, 0);
    for my $feat (@$feature) {
	my @part = split /\|\|/, $feat;
	if ($part[2] eq 'target_site_duplication') {
	    if ($tsd_ct > 0) {
		$five_pr_tsd = $part[4] - $part[3] + 1;
	    }
	    else {
		$three_pr_tsd = $part[4] - $part[3] + 1;
		$tsd_ct++;
	    }
	}
	elsif ($part[2] eq 'LTR_retrotransposon') {
	    ($ltr_sim) = ($part[8] =~ /ltr_similarity=?\s+?\"?(\d+.\d+)\"?/); 
	}
	$has_pbs = 1 if $part[2] eq 'primer_binding_site';
	$has_ppt = 1 if $part[2] eq 'RR_tract';
	$has_ir  = 1 if $part[2] eq 'inverted_repeat';
	$has_pdoms++ if $part[2] eq 'protein_match';
    }

    die "\nERROR: 'ltr_similarity' is not defined in GFF3. This is a bug, please report it.\n"
	unless defined $ltr_sim;
    $tsd_eq = 1 if $five_pr_tsd == $three_pr_tsd;

    my $ltr_score = sum($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq);
    return ($ltr_score, $ltr_sim);
}

sub _filter_compound_elements {
    my $self = shift;
    my ($features, $fasta) = @_;
    my (@pdoms, %classII);
    my $is_gypsy   = 0;
    my $is_copia   = 0;
    my $has_tpase  = 0;
    my $has_pdoms  = 0;
    my $len_thresh = 25000; # are elements > 25kb real? probably not

    my ($allct, $curct) = (0, 0);
    my ($gyp_cop_filtered, $dup_pdoms_filtered, $len_filtered, 
	$n_perc_filtered, $classII_filtered) = (0, 0, 0, 0, 0);
    
    for my $source (keys %$features) {
	for my $ltr (keys %{$features->{$source}}) {
	    $allct++;
	}
    }
    
    for my $source (keys %$features) {
	for my $ltr (keys %{$features->{$source}}) {
	    $curct++;
	    my ($rreg, $s, $e, $l) = split /\./, $ltr;
	    
	    for my $feat (@{$features->{$source}{$ltr}}) {
		my @feats = split /\|\|/, $feat;
		$feats[8] =~ s/\s\;\s/\;/g;
		$feats[8] =~ s/\s+/=/g;
		$feats[8] =~ s/\s+$//;
		$feats[8] =~ s/=$//;
		$feats[8] =~ s/=\;/;/g;
		$feats[8] =~ s/\"//g;

		if ($feats[2] =~ /protein_match/) {
		    $has_pdoms = 1;
		    @pdoms = ($feats[8] =~ /name=(\S+)/g);
		    if ($feats[8] =~ /name=RVT_1|name=Chromo/i) {
			$is_gypsy  = 1;
		    }
		    elsif ($feats[8] =~ /name=RVT_2/i) {
			$is_copia  = 1;
		    }
		    elsif ($feats[8] =~ /transpos|mule|dde|tnp_/i) {
			$has_tpase = 1;
		    }
		}
	    }
	    
	    if ($is_gypsy && $is_copia) {
		if (exists $features->{$source}{$ltr}) {
		    delete $features->{$source}{$ltr};
		    $gyp_cop_filtered++;
		}
	    }

	    if ($self->remove_tnp_domains && $has_tpase) {
		if (exists $features->{$source}{$ltr}) {
		    delete $features->{$source}{$ltr};
		    $classII_filtered++;
		}
	    }
	    
	    my %uniq;
	    if ($self->remove_dup_domains && $has_pdoms) {
		for my $element (@pdoms) {
		    $element =~ s/\;.*//;
		    next if $element =~ /chromo/i; # we expect these elements to be duplicated
		    if (exists $features->{$source}{$ltr} && $uniq{$element}++) {
			delete $features->{$source}{$ltr};
			$dup_pdoms_filtered++;# if $uniq{$element}++;
		    }
		}
	    }
	    
	    if ($l >= $len_thresh) {
		if (exists $features->{$source}{$ltr}) {
		    delete $features->{$source}{$ltr};
		    $len_filtered++;
		}
	    }

	    @pdoms = ();
	    $is_gypsy  = 0;
	    $is_copia  = 0;
	    $has_pdoms = 0;
	    $has_tpase = 0;
	}
    }

    my %stats = ( compound_gyp_cop_filtered   => $gyp_cop_filtered, 
		  length_filtered             => $len_filtered );

    if ($self->remove_dup_domains) {
	$stats{dup_pdoms_filtered} = $dup_pdoms_filtered;
    }
    if ($self->remove_tnp_domains) {
	$stats{class_II_filtered}  = $classII_filtered;
    }
    
    return $features, \%stats;
}

sub _filterNpercent {
    my $self = shift;
    my ($source, $key, $fasta) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($samtools) = @{$config}{qw(samtools)};
    my $n_perc = 0;
    my ($element, $start, $end) = split /\|\|/, $key;
    my $tmp = $element.".fasta";    
    my $cmd = "$samtools faidx $fasta $source:$start-$end > $tmp";
    $self->run_cmd($cmd);
    
    my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
    while (my $seqobj = $seqio->next_seq) {
	my $seq = $seqobj->seq;
	my $ltrlength = $seqobj->length;
	die unless length($seq) == $ltrlength;
	my $n_count = ($seq =~ tr/Nn//);
	$n_perc  = sprintf("%.2f",$n_count/$ltrlength);
    }
    unlink $tmp;

    return $n_perc;
}

sub _get_source {
    my $self = shift;
    my ($ref) = @_;
    for my $feat (@$ref) {
	my @feats = split /\|\|/, $feat;
	return ($feats[0]);
    }
}

sub _index_ref {
    my $self = shift;
    my ($samtools, $fasta) = @_;
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

    perldoc Tephra::LTR::LTRRefine


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
