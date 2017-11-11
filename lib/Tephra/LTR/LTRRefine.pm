package Tephra::LTR::LTRRefine;

use 5.014;
use Moose;
use File::Spec;
use File::Find;
use File::Basename;
use File::Copy          qw(move);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use Cwd                 qw(abs_path);
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Set::IntervalTree;
use Path::Class::File;
use Carp 'croak';
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util', 
     'Tephra::Role::Logger';

=head1 NAME

Tephra::LTR::LTRRefine - Refine LTR predictions based on multiple criteria

=head1 VERSION

Version 0.09.5

=cut

our $VERSION = '0.09.5';
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

has logfile => (
      is        => 'ro',
      isa       => 'Str',
      predicate => 'has_logfile',
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

has domains_required => (
    is        => 'ro',
    isa       => 'Bool',
    predicate => 'has_domains_required',
    lazy      => 1,
    default   => 0,
);

has is_trim => (
    is        => 'ro',
    isa       => 'Bool',
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

    my (%features, %intervals, %coord_map);
    my ($source, $seq_id, $start, $end, $length, $region);

    open my $gffio, '<', $gff or die "\nERROR: Could not open file: $gff\n";

    while (my $line = <$gffio>) {
	chomp $line;
	next if $line =~ /^#/;
	my $feature = gff3_parse_feature( $line );
	if ($feature->{type} eq 'repeat_region') {
	    $region = @{$feature->{attributes}{ID}}[0];
	    ($seq_id, $start, $end)  = @{$feature}{qw(seq_id start end)};
	    $length = $end - $start + 1; 
	    my $region_key = join "||", $region."_".$pid_thresh, $start, $end, $length;
	    $features{$seq_id}{$region_key} = [];
	    #$coord_map{$region."_".$pid_thresh} = join "||", $seq_id, $start, $end;
	    $intervals{$region."_".$pid_thresh} = join "||", $start, $end, $length;
	}
	
	if ($feature->{type} ne 'repeat_region') {
	    my $id = @{$feature->{attributes}{ID}}[0];
	    my $parent = @{$feature->{attributes}{Parent}}[0];
	        
	    if ($feature->{start} >= $start && $feature->{end} <= $end) {
		my $region_key = join "||", $region."_".$pid_thresh, $start, $end, $length;
		if ($parent =~ /(repeat_region\d+)/) {
		    my $old_parent = $1;
		    my $new_parent = join "_", $old_parent, $pid_thresh;
		    $parent =~ s/$old_parent/$new_parent/;
		    $feature->{attributes}{Parent}[0] = $parent;
		}
		push @{$features{$seq_id}{$region_key}}, $feature;
	    }
	}
    }

    my ($filtered, $stats) = $self->filter_compound_elements(\%features);
    return ({ collected_features => $filtered, stats => $stats, intervals => \%intervals}); #, region_map => \%coord_map });
}

sub filter_compound_elements {
    my $self = shift;
    my ($features) = @_;
    my (@pdoms, %classII);
    my $is_gypsy   = 0;
    my $is_copia   = 0;
    my $has_tpase  = 0;
    my $has_pdoms  = 0;
    my $len_thresh = 25000; # are elements > 25kb real? probably not for most systems

    my ($allct, $curct) = (0, 0);
    my ($gyp_cop_filtered, $dup_pdoms_filtered, $len_filtered, 
	$n_perc_filtered, $classII_filtered, $pdom_filtered) = (0, 0, 0, 0, 0, 0);
    
    for my $source (keys %$features) {
	for my $ltr (keys %{$features->{$source}}) {
	    $allct++;
	}
    }
    
    for my $source (keys %$features) {
	for my $ltr (keys %{$features->{$source}}) {
	    $curct++;
	    my ($rreg, $start, $end, $length) = split /\|\|/, $ltr;

	    for my $feature (@{$features->{$source}{$ltr}}) {
		if ($feature->{type} eq 'protein_match') {
		    my $pdom_name = $feature->{attributes}{name}[0];
		    push @pdoms, $pdom_name;
		    $has_pdoms = 1;

		    if ($pdom_name =~ /RVT_1|Chromo/i) {
			$is_gypsy  = 1;
		    }
		    elsif ($pdom_name =~ /RVT_2/i) {
			$is_copia  = 1;
		    }
		    elsif ($pdom_name =~ /transpos(?:ase)?|mule|(?:dbd|dde)?_tnp_(?:hat)?|duf4216/i) {
			# same regex when filtering elements for LTR family classification
			$has_tpase = 1;
		    }
		}
	    }
	            
	    if ($is_gypsy && $is_copia && exists $features->{$source}{$ltr}) {
		delete $features->{$source}{$ltr};
		$gyp_cop_filtered++;
	    }

	    if ($self->remove_tnp_domains && $has_tpase && exists $features->{$source}{$ltr}) {
		delete $features->{$source}{$ltr};
		$classII_filtered++;
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
	            
	    if ($length >= $len_thresh && exists $features->{$source}{$ltr}) {
		delete $features->{$source}{$ltr};
		$len_filtered++;
	    }

	    if ($self->domains_required && !$has_pdoms) {
		delete $features->{$source}{$ltr};
		$pdom_filtered++;
	    }

	    @pdoms = ();
	    $is_gypsy  = 0;
	    $is_copia  = 0;
	    $has_pdoms = 0;
	    $has_tpase = 0;
	}
    }

    my %stats = ( compound_gyp_cop_filtered => $gyp_cop_filtered, 
		  length_filtered           => $len_filtered );

    if ($self->remove_dup_domains) {
	$stats{dup_pdoms_filtered} = $dup_pdoms_filtered;
    }
    if ($self->remove_tnp_domains) {
	$stats{class_II_filtered} = $classII_filtered;
    }
    if ($self->domains_required) {
        $stats{coding_domain_filtered} = $pdom_filtered;
    }
    
    return ($features, \%stats);
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
	    my ($reg, $start, $end, $length) = split /\|\|/, $rregion;
	    $tree->insert($reg, $start, $end);
	    $chr_intervals{$source} = $tree;
	}
    }

    for my $source (keys %$partfeatures) {
	for my $rregion (keys %{$partfeatures->{$source}}) {
	    my (%scores, %sims);

	    if (exists $chr_intervals{$source}) {
		my ($reg, $start, $end, $length) = split /\|\|/, $rregion;
		my $res = $chr_intervals{$source}->fetch($start, $end);

		next unless defined $partfeatures->{$source}{$rregion};
		## This is the Tephra algorithm for scoring the quality of LTR-RT elements.
		##
		## The criteria are:
		## 1) the number of protein domains, 
		## 2) presence of TSDs,
		## 3) equal length TSDs, 
		## 4) polypurine (RR) tract,
		## 5) inverted repeats,
		## 6) primer binding site.
		##
		## If two elements are equal in these respects we take the element with the highest
		## LTR similarity.
		my ($score99, $sim99) = $self->summarize_features($partfeatures->{$source}{$rregion});
		
		$scores{$source}{$rregion} = $score99;
		$sims{$source}{$rregion} = $sim99;
		
		## Collect all features and LTR similary for overlapping elements.
		for my $overlap (@$res) {
		    my ($s, $e, $l) = split /\|\|/, $intervals->{$overlap};
		    my $region_key = join "||", $overlap, $s, $e, $l;
		    next unless defined $allfeatures->{$source}{$region_key};
		    my ($score85, $sim85) = $self->summarize_features($allfeatures->{$source}{$region_key});

		    $scores{$source}{$region_key} = $score85;
		    $sims{$source}{$region_key} = $sim85;
		}

		## Compute the best element based on the scores.
		if (@$res > 0) {
		    my $best_element 
			= $self->get_ltr_score_dups(\%scores, \%sims, $allfeatures, $partfeatures);
		    push @best_elements, $best_element;
		}
	    }
	}
    }
    
    return \@best_elements;
}

sub summarize_features {
    my $self = shift;
    my ($feature) = @_;
    my ($three_pr_tsd, $five_pr_tsd, $ltr_sim);
    my ($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq, $tsd_ct) = (0, 0, 0, 0, 0, 0);

    for my $feat (@$feature) {
	if ($feat->{type} eq 'target_site_duplication') {
	    if ($tsd_ct > 0) {
		$five_pr_tsd = $feat->{end} - $feat->{start} + 1;
	    }
	    else {
		$three_pr_tsd = $feat->{end} - $feat->{start} + 1;
		$tsd_ct++;
	    }
	}
	elsif ($feat->{type} eq 'LTR_retrotransposon') {
	    $ltr_sim = $feat->{attributes}{ltr_similarity}[0];
	}
	$has_pbs = 1 if $feat->{type} eq 'primer_binding_site';
	$has_ppt = 1 if $feat->{type} eq 'RR_tract';
	$has_ir  = 1 if $feat->{type} eq 'inverted_repeat';
	$has_pdoms++ if $feat->{type} eq 'protein_match';
    }

    croak "\nERROR: 'ltr_similarity' is not defined in GFF3. This is a bug, please report it.\n"
	unless defined $ltr_sim;
    $tsd_eq = 1 if $five_pr_tsd == $three_pr_tsd;

    my $ltr_score = sum($has_pbs, $has_ppt, $has_pdoms, $has_ir, $tsd_eq);
    return ($ltr_score, $ltr_sim);
}

sub get_ltr_score_dups {
    my $self = shift;
    my ($scores, $sims, $allfeatures, $partfeatures) = @_;
    my %sccounts;
    my %sicounts;
    my %best_element;
    my ($score_best, $sims_best) = (0, 0);

    ## These structures are to find out if we have multiple elements with
    ## the same features.
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

    my $chr = (keys %$scores)[0];
    my ($best_score_key, $best_sim_key, $best_sim, $best_score);
    for my $src (keys %$scores) {
	for my $rregion (keys %{$scores->{$src}}) {
	    if (defined $best_score) {
		$best_score_key = $rregion if $scores->{$src}{$rregion} >= $best_score;
	    }
	    else {
		$best_score = $scores->{$src}{$rregion};
		$best_score_key = $rregion;
	    }
	}
    }

    for my $src (keys %$sims) {
	for my $rregion (keys %{$sims->{$src}}) {
	    if (defined $best_sim) {
		$best_sim_key = $rregion if $sims->{$src}{$rregion} >= $best_sim;
	    }
	    else {
		$best_sim = $sims->{$src}{$rregion};
		$best_sim_key = $rregion;
	    }
	}
    }

    ## First, evaluate whether there is only one best element by score and similarity.
    ## Then, if there are multiple elements with the same best score/similarity, pick one
    ## based on the similarity.
    for my $source (keys %$scores) {
	if (@{$sccounts{ $scores->{$source}{$best_score_key} }} == 1 &&
	    @{$sicounts{ $sims->{$source}{$best_sim_key} }} == 1  &&
	    $best_score_key eq $best_sim_key) {
	    $score_best = 1;
	    if (exists $partfeatures->{$source}{$best_score_key}) {
		$best_element{$source}{$best_score_key} = $partfeatures->{$source}{$best_score_key};
	    }
	    elsif (exists $allfeatures->{$source}{$best_score_key}) {
		$best_element{$source}{$best_score_key} = $allfeatures->{$source}{$best_score_key};
	    }
	    else {
		croak "\nERROR: Something went wrong....'$best_score_key' not found in hash. This is a bug, please report it.";
	    }
	}
	elsif (@{$sccounts{ $scores->{$source}{$best_score_key} }} >= 1 && 
	       @{$sicounts{ $sims->{$source}{$best_sim_key} }} >= 1) {
	    $sims_best = 1;
	    my $best   = @{$sicounts{ $sims->{$source}{$best_sim_key} }}[0];
	    if (exists $partfeatures->{$source}{$best}) {
		$best_element{$source}{$best} = $partfeatures->{$source}{$best};
	    }
	    elsif (exists $allfeatures->{$source}{$best}) {
		$best_element{$source}{$best} = $allfeatures->{$source}{$best};
	    }
	    else {
		croak "\nERROR: Something went wrong....'$best' not found in hash. This is a bug, please report it.";
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

sub reduce_features {
    my $self = shift;
    my ($feature_ref) = @_;

    my $fasta = $self->genome->absolute->resolve;
    my $index = $self->index_ref($fasta);

    my ($logfile, $log);
    if (defined $self->logfile) {
	$logfile = $self->logfile;
	$log = $self->get_tephra_logger($logfile);
    }
    else {
	my ($name, $path, $suffix) = fileparse($fasta, qr/\.[^.]*/);
	my $lname = $self->is_trim ? 'tephra_findtrims.log' : 'tephra_findltrs.log';
	$logfile = File::Spec->catfile( abs_path($path), $name.'_'.$lname );
	$log = $self->get_tephra_logger($logfile);
	say STDERR "\nWARNING: '--logfile' option not given so results will be appended to: $logfile.";
    }

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
	    my ($region, $start, $end, $length) = split /\|\|/, $element;
	    my $key = join "||", $region, $start, $end;
	    
	    my $n_perc = $self->_filterNpercent($source, $key, $index);
	    if ($n_perc >= $n_thresh) {
		#say STDERR "=====> Over threshold: $n_perc";
		delete $best_features{$source}{$element};
		$n_perc_filtered++;
		$comb--;
	    }
	}
    }

    $best_stats{n_perc_filtered} = $n_perc_filtered;

    for my $s (keys %best_stats) {
	my $l = 40 - length($s);
	my $pad = ' ' x $l;
	$log->info("Results - Number of elements filtered by '$s':$pad",$best_stats{$s});
    }

    $log->info("Results - Number of elements found with 'relaxed' constraints:                       $all");
    $log->info("Results - Number of elements found with 'strict' constraints:                        $part");
    $log->info("Results - Number of 'best' elements that were overlapping in these two data sets:    $best");
    $log->info("Results - Number of 'combined' non-overlapping elements:                             $comb");

    return \%best_features;
}

sub sort_features {
    my $self = shift;
    my ($feature_ref) = @_;

    my ($gff, $combined_features) = @{$feature_ref}{qw(gff combined_features)};
    my $fasta = $self->genome->absolute->resolve;
    my $index = $self->index_ref($fasta);
    #my $logfile = $self->logfile;
    #my $log     = $self->get_logger($logfile);

    my ($outfile, $outfasta);
    if ($self->has_outfile) {
	$outfile = $self->outfile;
	my ($name, $path, $suffix) = fileparse($outfile, qr/\.[^.]*/);
	$outfasta = File::Spec->catfile( abs_path($path), $name.'.fasta' );
    }
    else {
	my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
	my $label = $self->is_trim ? 'tephra_trims' : 'tephra_ltrs';
	$name =~ s/ltrdigest85/$label/;
	$outfasta = File::Spec->catfile( abs_path($path), $name.'_combined_filtered.fasta' );
	$outfile  = File::Spec->catfile( abs_path($path), $name.'_combined_filtered.gff3' );
    }

    my ($logfile, $log);
    if (defined $self->logfile) {
        $logfile = $self->logfile;
        $log = $self->get_tephra_logger($logfile);
    }
    else {
        my ($name, $path, $suffix) = fileparse($fasta, qr/\.[^.]*/);
        my $lname = $self->is_trim ? 'tephra_findtrims.log' : 'tephra_findltrs.log';
        $logfile = File::Spec->catfile( abs_path($path), $name.'_'.$lname );
        $log = $self->get_tephra_logger($logfile);
	#no need to print warning twice
        #say STDERR "\nWARNING: '--logfile' option not given so results will be appended to: $logfile.";
    }

    open my $ofas, '>>', $outfasta or die "\nERROR: Could not open file: $outfasta\n";

    my ($elem_tot, $count) = (0, 1);
    if (defined $combined_features) {
	#dd $combined_features;
	open my $ogff, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

	my ($header, %features);
	open my $in, '<', $gff or die "\nERROR: Could not open file: $gff\n";
	while (my $line = <$in>) {
	    chomp $line;
	    if ($line =~ /^##\w+/) {
		$header .= $line."\n";
	    }
	    else {
		last;
	    }
	}
	close $in;

	if (not defined $header) {
	    say STDERR "\nWARNING: Could not get sequence region from $gff.\n";
	}
	else {
	    chomp $header;
	    say $ogff $header;
	}

	for my $chromosome (nsort keys %$combined_features) {
	    for my $ltr (nsort_by { m/repeat_region\d+\_\d+\|\|(\d+)\|\|\d+\|\|\d+/ and $1 }
			 keys %{$combined_features->{$chromosome}}) {
		my ($rreg, $rreg_start, $rreg_end, $rreg_length) = split /\|\|/, $ltr;
		my $new_rreg = $rreg;
		$new_rreg =~ s/\d+.*/$count/;

		my ($first) = @{$combined_features->{$chromosome}{$ltr}}[0];
		my ($source, $strand) = @{$first}{qw(source strand)};

		say $ogff join "\t", $chromosome, $source, 'repeat_region', 
		    $rreg_start, $rreg_end, '.', $strand, '.', "ID=$new_rreg";

		my ($parent, $gff_feats);
		for my $entry (@{$combined_features->{$chromosome}{$ltr}}) {
		    if ($entry->{type} eq 'LTR_retrotransposon') {
			$elem_tot++;

			my ($start, $end) = @{$entry}{qw(start end)};
			my $elem = $entry->{attributes}{ID}[0];
			$elem =~ s/\d+.*/$count/;

			my $id = join "_", $elem, $chromosome, $start, $end;
			#$id =~ s/LTR_/RLT_TRIM_/g
			$id =~ s/LTR_/TRIM_/g
			    if $self->is_trim;

			$self->_get_ltr_range($index, $id, $chromosome, $start, $end, $ofas);
			$entry->{attributes}{ID}[0] = $elem;
			$entry->{attributes}{Parent}[0] = $new_rreg;
			$parent = $elem;
		    }
		    elsif ($entry->{type} eq 'target_site_duplication' || $entry->{type} eq 'inverted_repeat') {
			$entry->{attributes}{Parent}[0] = $new_rreg;
		    }
		    else {
			$entry->{attributes}{Parent}[0] = $parent;
		    }
		    my $gff3_str = gff3_format_feature($entry);
		    $gff_feats .= $gff3_str;
		}
		$gff_feats =~ s/LTR_/TRIM_/g 
		    if $self->is_trim;

		print $ogff $gff_feats;
		$count++;
		undef $gff_feats;
	    }
	}
	close $ogff;

	#"Number of 'combined' non-overlapping elements:                              439"
	#say STDERR "\nTotal elements written: $elem_tot";
	my $pad = ' ' x 51;	       
	$log->info("Results - Total elements written:$pad",$elem_tot);
    }
    else {
	open my $in, '<', $gff, or die "\nERROR: Could not open file: $gff\n";
	while (my $line = <$in>) {
	    chomp $line;
	    next if $line =~ /^#/;
	    my $feature = gff3_parse_feature( $line );

	    if ($feature->{type} eq 'LTR_retrotransposon') {
		$elem_tot++;
		my ($chromosome, $start, $end) = @{$feature}{qw(seq_id start end)};

		my $elem = $feature->{attributes}{ID}[0];
		$elem =~ s/\d+.*//;
		$elem .= $count;
		#$elem =~ s/LTR_/RLT_TRIM_/g 
		$elem =~ s/LTR_/TRIM_/g
		    if $self->is_trim;

		my $id = join "_", $elem, $chromosome, $start, $end;
		$self->_get_ltr_range($index, $id, $chromosome, $start, $end, $ofas);
	    }
	}
	close $in;

	move $gff, $outfile or die "Move failed: $!";
	#say STDERR "\nTotal elements written: $elem_tot";
	#$log->info("Results - Total elements written: $elem_tot");
	my $pad = ' ' x 51;
        $log->info("Results - Total elements written:$pad",$elem_tot);
    }
    close $ofas;

    #say STDERR "\nTotal elements written: $elem_tot";
}

sub _get_ltr_range {
    my $self = shift;
    my ($index, $id, $chromosome, $start, $end, $ofh) = @_;
    my $location = "$chromosome:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);

    chomp $seq;
    $seq =~ s/^\s+|\s+$//g;
    $seq =~ s/.{60}\K/\n/g;
    say $ofh join "\n", ">".$id, $seq;
}

sub _filterNpercent {
    my $self = shift;
    my ($source, $key, $index) = @_;

    my ($element, $start, $end) = split /\|\|/, $key;
    my $location = "$source:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);

    my $n_count = ($seq =~ tr/Nn//);
    my $n_perc  = sprintf("%.2f",$n_count/$length);

    return $n_perc;
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

    perldoc Tephra::LTR::LTRRefine


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
