 #!/usr/bin/env perl

use 5.020;
use warnings;
use autodie;
use File::Basename;
use Statistics::Descriptive;
use Bio::Tools::GFF;
use Sort::Naturally     qw(nsort);
use List::Util          qw(sum max);
use List::MoreUtils     qw(first_index);
use List::UtilsBy       qw(nsort_by);
use IPC::System::Simple qw(system);
use Set::IntervalTree;
use Data::Dump;
use Data::Printer;
use Try::Tiny;
use Getopt::Long;
use experimental 'signatures';

my $usage = "$0 -fg ltrdigest85.gff3 -pg ltrdigest99.gff3 -f fasta -og out.gff3";
my $fullgff;
my $partgff;
my $fasta;  
my $outfile;

GetOptions(
    'fg|fullgff85=s' => \$fullgff,
    'pg|partgff99=s' => \$partgff,
    'f|fasta=s'      => \$fasta,
    'og|outgff=s'    => \$outfile,
    );

if (!$fullgff || !$partgff || !$fasta || !$outfile) {
    say $usage;
    exit(1);
}

my ($all_feats, $all_stats, $intervals)  = collect_features($fullgff, $fasta, '85');
my ($part_feats, $part_stats, $part_int) = collect_features($partgff, $fasta, '99');
#dd $all_feats and exit;

my $best_elements = get_overlaps($all_feats, $part_feats, $intervals);

my $combined_features = reduce_features($all_feats, $part_feats, $best_elements, 
					$all_stats, $part_stats, $fasta);
sort_features($fullgff, $fasta, $combined_features, $outfile);

exit;
#
# methods
#
sub reduce_features ($all_feats, $part_feats, $best_elements, $all_stats, $part_stats, $fasta) {
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
    my $n_thresh = 0.30;

    for my $source (keys %best_features) {
	for my $element (keys %{$best_features{$source}}) {
	    $comb++;
	    my ($region, $start, $end, $length) = split /\./, $element;
	    my $key = join "||", $region, $start, $end;

	    my $n_perc = filterNpercent($source, $key, $fasta);
	    if ($n_perc >= $n_thresh) {
		#say STDERR "=====> Over thresh: $n_perc";
		delete $best_features{$source}{$element};
		$n_perc_filtered++;
		$comb--;
	    }
	}
    }

    say STDERR "all stats: ";
    p $all_stats;
    say STDERR "part stats: ";
    p $part_stats;
    say STDERR "best stats: ";
    p %best_stats;
    $best_stats{n_perc_filtered} = $n_perc_filtered;
    say STDERR "best stats: ";
    p %best_stats;
    

    say STDERR join q{ }, "filtered_type", "num_filtered";
    for my $s (keys %best_stats) {
	say STDERR join q{ }, $s, $best_stats{$s};
    }

    say STDERR join q{ }, "All", "part", "best", "combined";
    say STDERR join q{ }, $all, $part, $best, $comb;

    return \%best_features;
}

sub sort_features ($gff, $fasta, $combined_features, $outfile) {
    my ($name, $path, $suffix) = fileparse($outfile, qr/\.[^.]*/);
    my $outfasta = $name.".fasta";
    open my $ogff, '>', $outfile;
    open my $ofas, '>>', $outfasta;

    my ($header, %features);
    open my $in, '<', $gff;
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

    my $elem_tot = 0;
    for my $chromosome (nsort keys %$combined_features) {
	for my $ltr (nsort_by { m/repeat_region\d+\_\d+\.(\d+)\.\d+/ and $1 }
		     keys %{$combined_features->{$chromosome}}) {
	    my ($rreg, $rreg_start, $rreg_end, $rreg_length) = split /\./, $ltr;
	    my ($first) = @{$combined_features->{$chromosome}{$ltr}}[0];
	    my ($source, $strand) = (split /\|\|/, $first)[1,6];
	    say $ogff join "\t", $chromosome, $source, 'repeat_region', 
	        $rreg_start, $rreg_end, '.', $strand, '.', "ID=$rreg";
	    for my $entry (@{$combined_features->{$chromosome}{$ltr}}) {
		my @feats = split /\|\|/, $entry;
		$feats[8] =~ s/\s\;\s/\;/g;
		$feats[8] =~ s/\s+$//;
		$feats[8] =~ s/\"//g;
		$feats[8] =~ s/(\;\w+)\s/$1=/g;
		$feats[8] =~ s/\s;/;/;
		$feats[8] =~ s/^(\w+)\s/$1=/;
		say $ogff join "\t", @feats;

		if ($feats[2] eq 'LTR_retrotransposon') {
		    $elem_tot++;
		    my ($start, $end) = @feats[3..4];
		    my ($elem) = ($feats[8] =~ /(LTR_retrotransposon\d+)/);
		    my $id = $elem."_".$chromosome."_".$start."_".$end;
		    my $tmp = $elem.".fasta";
		    my $cmd = "samtools faidx $fasta $chromosome:$start-$end > $tmp";
		    try {
			system([0..5], $cmd);
		    }
		    catch {
			die "\nERROR: $cmd failed. Here is the exception: $_\n";
		    };
		    my @aux = undef;
		    my ($name, $comm, $seq, $qual);
		    open my $in, '<', $tmp;
		    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
			$seq =~ s/.{60}\K/\n/g;
			say $ofas join "\n", ">".$id, $seq;
		    }   
		    close $in;
		    unlink $tmp;
		} 
	    }
	}
    }
    close $ogff;
    close $ofas;

    say STDERR "Total elements written: $elem_tot";
}
    
sub collect_features ($gff, $fasta, $which) {
    my %intervals;
    
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
		$intervals{$region."_".$which} = join ".", $start, $end, $length;
		my $region_key = join ".", $region."_".$which, $start, $end, $length;
		my @feats = split /\t/, $feature->gff_string;
		if ($feats[8] =~ /(repeat_region\d+)/) {
		    my $old_parent = $1;
		    my $new_parent = $old_parent."_".$which;
		    $feats[8] =~ s/$old_parent/$new_parent/g;
		}
		push @{$features{$source}{$region_key}}, join "||", @feats; #split /\t/, $feature->gff_string;
	    }
	}
    }

    my ($filtered, $stats) = filter_compound_elements(\%features, $fasta);
    return $filtered, $stats, \%intervals;
}

sub get_overlaps ($allfeatures, $partfeatures, $intervals) {
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

		my ($score99, $sim99) = summarize_features($partfeatures->{$source}{$rregion});
		
		$scores{$source}{$rregion} = $score99;
		$sims{$source}{$rregion} = $sim99;
		
		## collect all features, then compare scores/sim...
		for my $over (@$res) {
		    my ($s, $e, $l) = split /\./, $intervals->{$over};
		    my $region_key = join ".", $over, $s, $e, $l;
		    my ($score85, $sim85) = summarize_features($allfeatures->{$source}{$region_key});
		    
		    $scores{$source}{$region_key} = $score85;
		    $sims{$source}{$region_key} = $sim85;
		}

		if (@$res > 0) {
		    my $best_element = get_ltr_score_dups(\%scores, \%sims, $allfeatures, $partfeatures);
		    push @best_elements, $best_element;
		}
	    }
	}
    }
    
    return \@best_elements;
}

sub get_ltr_score_dups ($scores, $sims, $allfeatures, $partfeatures) {
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
	    my $bscore = $scores->{$source}{$best_score_key};
	    my $bsim   = $sims->{$source}{$best_score_key};
	    if (exists $partfeatures->{$source}{$best_score_key}) {
		$best_element{$source}{$best_score_key} = $partfeatures->{$source}{$best_score_key};
	    }
	    elsif (exists $allfeatures->{$source}{$best_score_key}) {
		$best_element{$source}{$best_score_key} = $allfeatures->{$source}{$best_score_key};
	    }
	    else {
		say "\nERROR: Something went wrong....'$best_score_key' not found in hash. This is a bug. Exiting.";
		exit(1);
	    }
	}
	elsif (@{$sccounts{ $scores->{$source}{$best_score_key} }} >= 1 && 
	       @{$sicounts{ $sims->{$source}{$best_sim_key} }} >=  1) {
	    $sims_best = 1;
	    my $best = @{$sicounts{ $sims->{$source}{$best_sim_key} }}[0];
	    my $bscore = $scores->{$source}{$best_score_key};
	    my $bsim   = $sims->{$source}{$best_score_key};
	    if (exists $partfeatures->{$source}{$best}) {
		$best_element{$source}{$best_score_key} = $partfeatures->{$source}{$best};
	    }
	    elsif (exists $allfeatures->{$source}{$best}) {
		$best_element{$source}{$best_score_key} = $allfeatures->{$source}{$best};
	    }
	    else {
		say "\nERROR: Something went wrong....'$best' not found in hash. This is a bug. Exiting.";
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

sub summarize_features ($feature) {
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
	    ($ltr_sim) = ($part[8] =~ /ltr_similarity \"(\d+.\d+)\"/); 
	}
	$has_pbs = 1 if $part[2] eq 'primer_binding_site';
	$has_ppt = 1 if $part[2] eq 'RR_tract';
	$has_ir  = 1 if $part[2] eq 'inverted_repeat';
	$has_pdoms++ if $part[2] eq 'protein_match';
    }
    p $feature and exit(1) unless defined $ltr_sim;
    $tsd_eq = 1 if $five_pr_tsd == $three_pr_tsd;

    my $ltr_score = sum($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq);
    return ($ltr_score, $ltr_sim);
}

sub filter_compound_elements ($features, $fasta) {
    my @pdoms;
    my $is_gypsy   = 0;
    my $is_copia   = 0;
    my $has_tpase  = 0;
    my $has_pdoms  = 0;
    my $len_thresh = 25000; # are elements > 25kb real? probably not
    my $n_thresh   = 0.30;

    my ($allct, $curct) = (0, 0);
    my ($gyp_cop_filtered, $dup_pdoms_filtered, $len_filtered, 
	$n_perc_filtered, $classII_filtered) = (0, 0, 0, 0, 0);
    
    for my $source (keys %$features) {
	for my $ltr (keys %{$features->{$source}}) {
	    $allct++;
	}
    }
    
    for my $source (keys %$features) {
	#for my $ltr (nsort_by { m/repeat_region(\d+)/ and $1 } keys %{$features->{$source}}) {
	for my $ltr (keys %{$features->{$source}}) {
	    $curct++;
	    my ($rreg, $s, $e, $l) = split /\./, $ltr;
	    #my $region = @{$features->{$source}{$ltr}}[0];
	    
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
		    elsif ($feats[8] =~ /transposase/i) {
			$has_tpase = 1;
		    }
		}
	    }
	    
	    if ($is_gypsy && $is_copia) {
		delete $features->{$source}{$ltr};
		$gyp_cop_filtered++;
	    }

	    if ($has_tpase) {
		delete $features->{$source}{$ltr};
		$classII_filtered++;
	    }
	    
	    my %uniq;
	    if ($has_pdoms) {
		for my $element (@pdoms) {
		    $element =~ s/\;.*//;
		    next if $element =~ /chromo/i; # we expect these elements to be duplicated
		    delete $features->{$source}{$ltr} && $dup_pdoms_filtered++ if $uniq{$element}++;
		}
	    }
	    
	    if ($l >= $len_thresh) {
		delete $features->{$source}{$ltr};
		$len_filtered++;
	    }

	    @pdoms = ();
	    $is_gypsy  = 0;
	    $is_copia  = 0;
	    $has_pdoms = 0;
	    $has_tpase = 0;
	}
    }

    my %stats = ( gyp_cop_filtered   => $gyp_cop_filtered, 
		  len_filtered       => $len_filtered, 
		  dup_pdoms_filtered => $dup_pdoms_filtered,
	          class_II_filtered  => $classII_filtered );

    return $features, \%stats;
}

sub filterNpercent ($source, $key, $fasta) {
    my $n_perc = 0;
    my ($element, $start, $end) = split /\|\|/, $key;
    my $tmp = $element.".fasta";
    my $cmd = "samtools faidx $fasta $source:$start-$end > $tmp";
    try { 
	system([0..5], $cmd);
    }
    catch {
	die "\nERROR: $cmd failed. Here is the exception: $_\n";
    };
    
    my @aux = undef;
    my ($name, $comm, $seq, $qual);
    open my $in, '<', $tmp;
    while (($name, $comm, $seq, $qual) = readfq(\*$in, \@aux)) {
	my $ltrlength = length($seq);
	my $n_count = ($seq =~ tr/Nn//);
	$n_perc  = sprintf("%.2f",$n_count/$ltrlength);
	#say STDERR join q{ }, $source, $ltr, $start, $end, $ltrlength, $n_perc;
    }
    close $in;
    unlink $tmp;

    return $n_perc;
}

sub get_source ($ref) {
    for my $feat (@$ref) {
	my @feats = split /\|\|/, $feat;
	return ($feats[0]);
    }
}

sub readfq {
    my ($fh, $aux) = @_;
    @$aux = [undef, 0] if (!@$aux);
    return if ($aux->[1]);
    if (!defined($aux->[0])) {
        while (<$fh>) {
            chomp;
            if (substr($_, 0, 1) eq '>' || substr($_, 0, 1) eq '@') {
                $aux->[0] = $_;
                last;
            }
        }
        if (!defined($aux->[0])) {
            $aux->[1] = 1;
            return;
        }
    }
    my ($name, $comm);
    defined $_ && do {
        ($name, $comm) = /^.(\S+)(?:\s+)(\S+)/ ? ($1, $2) : 
	                 /^.(\S+)/ ? ($1, '') : ('', '');
    };
    my $seq = '';
    my $c;
    $aux->[0] = undef;
    while (<$fh>) {
        chomp;
        $c = substr($_, 0, 1);
        last if ($c eq '>' || $c eq '@' || $c eq '+');
        $seq .= $_;
    }
    $aux->[0] = $_;
    $aux->[1] = 1 if (!defined($aux->[0]));
    return ($name, $comm, $seq) if ($c ne '+');
    my $qual = '';
    while (<$fh>) {
        chomp;
        $qual .= $_;
        if (length($qual) >= length($seq)) {
            $aux->[0] = undef;
            return ($name, $comm, $seq, $qual);
        }
    }
    $aux->[1] = 1;
    return ($name, $seq);
}
