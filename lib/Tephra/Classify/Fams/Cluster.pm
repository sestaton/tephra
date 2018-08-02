package Tephra::Classify::Fams::Cluster;

use 5.014;
use Moose::Role;
use MooseX::Types::Path::Class;
use Sort::Naturally;
use Number::Range;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Bio::GFF3::LowLevel qw(gff3_parse_feature);
use List::Util          qw(min max);
use Time::HiRes         qw(gettimeofday);
use File::Path          qw(make_path);
use Cwd                 qw(abs_path);
use List::UtilsBy       qw(nsort_by);
use Parallel::ForkManager;
use Carp 'croak';
use Try::Tiny;
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::Any';

=head1 NAME

     Tephra::Classify::Fams::Cluster - Helper role with clustering routines for LTR/TIR classification

=head1 VERSION

Version 0.12.0

=cut

our $VERSION = '0.12.0';
$VERSION = eval $VERSION;

has debug => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

sub extract_ltr_features {
    my $self = shift;
    my $fasta = $self->genome->absolute->resolve;
    my $dir   = $self->outdir->absolute->resolve;
    my ($infile) = @_;
    
    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $type = ($name =~ /(?:gypsy|copia|unclassified)$/i);
    croak "\n[ERROR]: Unexpected input. Should match /(?:gypsy|copia|unclassified)$/i. Exiting."
        unless defined $type;

    my $resdir = File::Spec->catdir($dir, $name);
    unless ( -d $resdir ) {
        make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $comp = File::Spec->catfile($resdir, $name.'_complete.fasta');
    my $ppts = File::Spec->catfile($resdir, $name.'_ppt.fasta');
    my $pbs  = File::Spec->catfile($resdir, $name.'_pbs.fasta');
    my $five_pr_ltrs  = File::Spec->catfile($resdir, $name.'_5prime-ltrs.fasta');
    my $three_pr_ltrs = File::Spec->catfile($resdir, $name.'_3prime-ltrs.fasta');

    open my $allfh, '>>', $comp or die "\n[ERROR]: Could not open file: $comp\n";
    open my $pptfh, '>>', $ppts or die "\n[ERROR]: Could not open file: $ppts\n";
    open my $pbsfh, '>>', $pbs or die "\n[ERROR]: Could not open file: $pbs\n";
    open my $fivefh, '>>', $five_pr_ltrs or die "\n[ERROR]: Could not open file: $five_pr_ltrs\n";
    open my $threfh, '>>', $three_pr_ltrs or die "\n[ERROR]: Could not open file: $three_pr_ltrs\n";

    open my $gffio, '<', $infile or die "\n[ERROR]: Could not open file: $infile\n";

    my (%feature, %ltrs, %coord_map, %seen);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
        my $feature = gff3_parse_feature( $line );

        if ($feature->{type} =~ /(?:LTR|TRIM|LARD)_retrotransposon/) {
            my $elem_id = @{$feature->{attributes}{ID}}[0];
            my ($start, $end) = @{$feature}{qw(start end)};
            my $key = join "||", $elem_id, $start, $end;
            $ltrs{$key}{'full'} = join "||", @{$feature}{qw(seq_id type start end)};
            $coord_map{$elem_id} = join "||", @{$feature}{qw(seq_id start end)};
	}
	if ($feature->{type} eq 'long_terminal_repeat') {
            my $parent = @{$feature->{attributes}{Parent}}[0];
            my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
                my ($seq_id, $type, $start, $end, $strand) = 
                    @{$feature}{qw(seq_id type start end strand)};
                $strand //= '?';
                my $ltrkey = join "||", $seq_id, $type, $start, $end, $strand;
                push @{$ltrs{$pkey}{'ltrs'}}, $ltrkey unless exists $seen{$ltrkey};
                $seen{$ltrkey} = 1;
            }
        }
        elsif ($feature->{type} eq 'primer_binding_site') {
            my $name = @{$feature->{attributes}{trna}}[0];
            my $parent = @{$feature->{attributes}{Parent}}[0];
            my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
                $ltrs{$pkey}{'pbs'} =
                    join "||", @{$feature}{qw(seq_id type)}, $name, @{$feature}{qw(start end)};
            }
        }
        elsif ($feature->{type} eq 'protein_match') {
	    my $name = @{$feature->{attributes}{name}}[0];
            my $parent = @{$feature->{attributes}{Parent}}[0];
            my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
                my $pdomkey = join "||", @{$feature}{qw(seq_id type)}, $name, @{$feature}{qw(start end strand)};
                #push @{$ltrs{$pkey}{'pdoms'}{$name}}, $pdomkey unless exists $seen{$pdomkey};
		push @{$ltrs{$pkey}{'pdoms'}}, $pdomkey unless exists $seen{$pdomkey};
                $seen{$pdomkey} = 1;
            }
	}
	elsif ($feature->{type} eq 'RR_tract') {
            my $parent = @{$feature->{attributes}{Parent}}[0];
            my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
                $ltrs{$pkey}{'ppt'} =
                    join "||", @{$feature}{qw(seq_id type start end)};
            }
        }
    }
    close $gffio;
    #dd \%ltrs;

    my (%pdoms, %seen_pdoms, %lrange);
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
        my ($element, $rstart, $rend) = split /\|\|/, $ltr;
        # full element
        my ($source, $prim_tag, $fstart, $fend) = split /\|\|/, $ltrs{$ltr}{'full'};
	my $fullid = join "_", $element, $source, $fstart, $fend;
	$self->write_element_parts($index, $source, $fstart, $fend, $allfh, $fullid);

        # pbs
        if ($ltrs{$ltr}{'pbs'}) {
            my ($pbssource, $pbstag, $trna, $pbsstart, $pbsend) = split /\|\|/, $ltrs{$ltr}{'pbs'};
	    my $pbsid = join "_", $element, $pbssource, $pbsstart, $pbsend;
	    $self->write_element_parts($index, $pbssource, $pbsstart, $pbsend, $pbsfh, $pbsid);
        }
        # ppt
        if ($ltrs{$ltr}{'ppt'}) {
            my ($pptsource, $ppttag, $pptstart, $pptend) = split /\|\|/, $ltrs{$ltr}{'ppt'};
	    my $pptid = join "_", $element, $pptsource, $pptstart, $pptend;
	    $self->write_element_parts($index, $pptsource, $pptstart, $pptend, $pptfh, $pptid);
	}
	
	# ltrs
	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
            my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
	    my $ltrid = join "_", $element, $src, $s, $e;
            if ($ltrct) {
		$self->write_element_parts($index, $src, $s, $e, $fivefh, $ltrid);
                $ltrct = 0;
            }
            else {
		$self->write_element_parts($index, $src, $s, $e, $threfh, $ltrid);
                $ltrct++;
            }
        }

	# pdoms
	if ($ltrs{$ltr}{'pdoms'}) {
 	    for my $ltr_repeat (@{$ltrs{$ltr}{'pdoms'}}) {
		my ($src, $pdomtag, $name, $s, $e, $str) = split /\|\|/, $ltr_repeat;
		#"Ha10||protein_match||UBN2||132013916||132014240|+",
		next if $name =~ /transpos(?:ase)?|mule|(?:dbd|dde)?_tnp_(?:hat)?|duf4216/i; 
		# The above is so we do not classify elements based domains derived from, or belonging to, DNA transposons
		push @{$pdoms{$src}{$element}}, join "||", $name, $s, $e, $str;
		push @{$lrange{$src}{$element}{$name}}, "$s..$e";
            }
        }
    }
    close $allfh;
    close $pptfh;
    close $pbsfh;
    close $fivefh;
    close $threfh;

    my $pdom_fam_map = $self->merge_overlapping_hits($index, $resdir, \%pdoms, \%lrange);

    for my $file ($comp, $ppts, $pbs, $five_pr_ltrs, $three_pr_ltrs) {
        unlink $file if ! -s $file;
    }

    return ($resdir, $pdom_fam_map);
}

sub extract_tir_features {
    my $self = shift;
    my $fasta = $self->genome->absolute->resolve;
    my $dir   = $self->outdir->absolute->resolve;
    my ($infile) = @_;
    
    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $type = ($name =~ /(?:cacta|mariner|mutator|hat|unclassified|mite)$/i);
    croak "\n[ERROR]: Unexpected input. Should match /cacta|mariner|mutator|hat|unclassified|mite$/i. Exiting."
        unless defined $type;

    my $resdir = File::Spec->catdir($dir, $name);
    unless ( -d $resdir ) {
        make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $comp = File::Spec->catfile($resdir, $name.'_complete.fasta');
    my $five_pr_tirs  = File::Spec->catfile($resdir, $name.'_5prime-tirs.fasta');
    my $three_pr_tirs = File::Spec->catfile($resdir, $name.'_3prime-tirs.fasta');

    open my $allfh, '>>', $comp or die "\n[ERROR]: Could not open file: $comp\n";
    open my $fivefh, '>>', $five_pr_tirs or die "\n[ERROR]: Could not open file: $five_pr_tirs\n";
    open my $threfh, '>>', $three_pr_tirs or die "\n[ERROR]: Could not open file: $three_pr_tirs\n";

    open my $gffio, '<', $infile or die "\n[ERROR]: Could not open file: $infile\n";

    my (%feature, %tirs, %coord_map, %seen);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
        my $feature = gff3_parse_feature( $line );

        if ($feature->{type} =~ /terminal_inverted_repeat_element|MITE/) {
            my $elem_id = @{$feature->{attributes}{ID}}[0];
            my ($start, $end) = @{$feature}{qw(start end)};
            my $key = join "||", $elem_id, $start, $end;
            $tirs{$key}{'full'} = join "||", @{$feature}{qw(seq_id type start end)};
            $coord_map{$elem_id} = join "||", @{$feature}{qw(seq_id start end)};
        }
        if ($feature->{type} eq 'terminal_inverted_repeat') {
            my $parent = @{$feature->{attributes}{Parent}}[0];
            my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
                my ($seq_id, $type, $start, $end, $strand) = 
                    @{$feature}{qw(seq_id type start end strand)};
                $strand //= '?';
                my $tirkey = join "||", $seq_id, $type, $start, $end, $strand;
                push @{$tirs{$pkey}{'tirs'}}, $tirkey unless exists $seen{$tirkey};
                $seen{$tirkey} = 1;
            }
        }
        elsif ($feature->{type} eq 'protein_match') {
	    my $name = @{$feature->{attributes}{name}}[0];
            my $parent = @{$feature->{attributes}{Parent}}[0];
            my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
            if ($seq_id eq $feature->{seq_id}) {
                my $pdomkey = join "||", @{$feature}{qw(seq_id type)}, $name, @{$feature}{qw(start end strand)};
                #push @{$tirs{$pkey}{'pdoms'}{$name}}, $pdomkey unless exists $seen{$pdomkey};
		push @{$tirs{$pkey}{'pdoms'}}, $pdomkey unless exists $seen{$pdomkey};
                $seen{$pdomkey} = 1;
            }
        }
    }
    close $gffio;

    my (%pdoms, %seen_pdoms, %lrange);
    my $tirct = 0;
    for my $tir (sort keys %tirs) {
        my ($element, $rstart, $rend) = split /\|\|/, $tir;
        # full element
        my ($source, $prim_tag, $fstart, $fend) = split /\|\|/, $tirs{$tir}{'full'};
	my $fullid = join "_", $element, $source, $fstart, $fend;
	$self->write_element_parts($index, $source, $fstart, $fend, $allfh, $fullid);

	# tirs
	my $partid;
        for my $tir_repeat (@{$tirs{$tir}{'tirs'}}) {
            my ($src, $tirtag, $s, $e, $strand) = split /\|\|/, $tir_repeat;
	    my $partid = join "_", $element, $src, $s, $e;

            if ($tirct) {
		$self->write_element_parts($index, $src, $s, $e, $fivefh, $partid);
                $tirct = 0;
            }
            else {
		$self->write_element_parts($index, $src, $s, $e, $threfh, $partid);
                $tirct++;
            }
        }

	# pdoms
	if ($tirs{$tir}{'pdoms'}) {
	    for my $tir_repeat (@{$tirs{$tir}{'pdoms'}}) {
		my ($src, $pdomtag, $name, $s, $e, $str) = split /\|\|/, $tir_repeat;
		#"Ha10||protein_match||UBN2||132013916||132014240|+",
		next unless $name =~ /transpos(?:ase)?|mule|(?:dbd|dde)?_tnp_(?:hat)?|duf4216/i; 
		#next if $model_name =~ /rve|rvt|rvp|gag|chromo|rnase|athila|zf/i)
		# The above is so we do not classify elements based domains derived from or belonging to DNA transposons
		push @{$pdoms{$src}{$element}}, join "||", $name, $s, $e, $str;
		push @{$lrange{$src}{$element}{$name}}, "$s..$e";
            }
        }
    }
    close $allfh;
    close $fivefh;
    close $threfh;
    
    my $pdom_fam_map = $self->merge_overlapping_hits($index, $resdir, \%pdoms, \%lrange);

    for my $file ($comp, $five_pr_tirs, $three_pr_tirs) {
        unlink $file if ! -s $file;
    }

    return ($resdir, $pdom_fam_map);
}

sub merge_overlapping_hits {
    my $self = shift;
    my ($index, $resdir, $pdoms, $lrange) = @_;

    #dd $pdoms and exit;
    ## This is where we merge overlapping hits in a chain and concatenate non-overlapping hits
    ## to create a single domain sequence for each element
    my (%pdom_fam_map, %element_map, %seen, %doms);
    for my $src (keys %$pdoms) {
        for my $element (keys %{$pdoms->{$src}}) {
            my ($pdom_name, $pdom_s, $pdom_e, $str);
	    
	    for my $pdom_type (@{$pdoms->{$src}{$element}}) {
		my (%seqs, $union);                                                                                               
		my ($pdom_name, $pdom_start, $pdom_eend, $strand) = split /\|\|/, $pdom_type;

                my $pdom_file = File::Spec->catfile( abs_path($resdir), $pdom_name.'_pdom.fasta' );        
                open my $fh, '>>', $pdom_file or die "\n[ERROR]: Could not open file: $pdom_file\n";
	    
		if (@{$lrange->{$src}{$element}{$pdom_name}} > 1) {
		    unless (exists $seen{$src}{$element}{$pdom_name}) {
			{
			    no warnings; # Number::Range warns on EVERY single interger that overlaps
			    my $range = Number::Range->new(@{$lrange->{$src}{$element}{$pdom_name}});
			    $union = $range->range;
			}

			for my $dom (split /\,/, $union) {
			    my ($ustart, $uend) = split /\.\./, $dom;
			    push @{$doms{$element}}, join "||", $pdom_name, $ustart, $uend;
			    #push @{$doms{$element}}, join "||", $ustart, $uend;
			    my ($seq, $length) = $self->get_full_seq($index, $src, $ustart, $uend);
			    my $k = join "||", $ustart, $uend;
			    $seqs{$k} = $seq;
			}
	
			$self->concat_pdoms($index, $src, $pdom_name, $element, \%seqs, $fh);
		    }
		    $seen{$src}{$element}{$pdom_name} = 1;
		}
		else {
		    my ($nuname, $nustart, $nuend, $str) = split /\|\|/, $pdom_type;
		    push @{$doms{$element}}, join "||", $nuname, $nustart, $nuend;
		    #push @{$doms{$element}}, join "||", $nustart, $nuend;
		    my $id = join "_", $element, $src, $nustart, $nuend;
		    $self->write_element_parts($index, $src, $nustart, $nuend, $fh, $id);
		}
		close $fh;
		%seqs = ();
		unlink $pdom_file if ! -s $pdom_file;
	    }
	    if (@{$doms{$element}}) {
		my @dom_order;
		for my $dom (nsort_by { m/\S+\|\|(\d+)\|\|\d+/ and $1 } @{$doms{$element}}) {
		    my ($dom_name, $start, $end) = split /\|\|/, $dom; 
		    push @dom_order, $dom_name."{$start-$end}";
		}
		$pdom_fam_map{$element}{pdoms} = join ",", @dom_order;
	    }
	    %doms = ();
	}
    }

    #dd \%pdom_fam_map; # and exit;
    return \%pdom_fam_map;
}

sub concat_pdoms {
    my $self = shift;
    my ($index, $src, $pdom_name, $element, $seqs, $fh_out) = @_;
    my $id = join "_", $element, $src."_"; #, $pdom_name."_";

    my $concat_seq;
    for my $dom (keys %$seqs) {
	my ($start, $end) = split /\|\|/, $dom;
	$id .= join "_", $start, $end."_";
	$concat_seq .= $seqs->{$dom};
    }
    $id =~ s/\_$//;
    $concat_seq =~ s/.{60}\K/\n/g;
    say $fh_out join "\n", ">$id", $concat_seq;

    return;
}

sub collect_feature_args {
    my $self = shift;
    my ($dir) = @_;
    my $tetype = $self->type;

    my (@fiveltrs, @threeltrs, @fivetirs, @threetirs, @ppt, @pbs, @pdoms, %vmatch_args);
    if ($tetype eq 'LTR') {
	find( sub { push @fiveltrs, $File::Find::name if -f and /5prime-ltrs.fasta$/ }, $dir);
	find( sub { push @threeltrs, $File::Find::name if -f and /3prime-ltrs.fasta$/ }, $dir);
	find( sub { push @ppt, $File::Find::name if -f and /ppts.fasta$/ }, $dir);
	find( sub { push @pbs, $File::Find::name if -f and /pbs.fasta$/ }, $dir);
	find( sub { push @pdoms, $File::Find::name if -f and /pdom.fasta$/ }, $dir);
	
	# ltr
	my $ltr5name = File::Spec->catfile( abs_path($dir), 'dbcluster-5primeseqs' );
	my $fiveargs = "-qspeedup 2 -dbcluster 80 0 $ltr5name -p -d -seedlength 30 ";
	$fiveargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 10000 -identity 80";
	$vmatch_args{fiveltr} = { seqs => \@fiveltrs, args => $fiveargs };

	my $ltr3name  = File::Spec->catfile( abs_path($dir), 'dbcluster-3primeseqs' );
	my $threeargs = "-qspeedup 2 -dbcluster 80 0 $ltr3name -p -d -seedlength 30 ";
	$threeargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 10000 -identity 80";
	$vmatch_args{threeltr} = { seqs => \@threeltrs, args => $threeargs };
	
	# pbs/ppt
	my $pbsname = File::Spec->catfile( abs_path($dir), 'dbcluster-pbs' );
	my $pbsargs = "-dbcluster 90 90 $pbsname -p -d -seedlength 5 -exdrop 2 ";
	$pbsargs .= "-l 3 -showdesc 0 -sort ld -best 10000 -identity 90";
	$vmatch_args{pbs} = { seqs => \@pbs, args => $pbsargs, prefixlen => 1 };
	
	my $pptname = File::Spec->catfile( abs_path($dir), 'dbcluster-ppt' );
	my $pptargs = "-dbcluster 90 90 $pptname -p -d -seedlength 5 -exdrop 2 ";
	$pptargs .= "-l 3 -showdesc 0 -sort ld -best 10000 -identity 90";
	$vmatch_args{ppt} = { seqs => \@ppt, args => $pptargs, prefixlen => 5 };
	
	# pdoms
	my $pdomname = File::Spec->catfile( abs_path($dir), 'dbcluster-pdoms' );
	my $pdomargs = "-qspeedup 2 -dbcluster 80 0 $pdomname -p -d -seedlength 30 -exdrop 3 ";
	$pdomargs .= "-l 40 -showdesc 0 -sort ld -best 10000";
	$vmatch_args{pdoms} = { seqs => \@pdoms, args => $pdomargs };
    }
    
    if ($tetype eq 'TIR') {
	find( sub { push @fivetirs, $File::Find::name if -f and /5prime-tirs.fasta$/ }, $dir);
	find( sub { push @threetirs, $File::Find::name if -f and /3prime-tirs.fasta$/ }, $dir);
	find( sub { push @pdoms, $File::Find::name if -f and /pdom.fasta$/ }, $dir);
	
	# tir
	my $tir5name = File::Spec->catfile( abs_path($dir), 'dbcluster-5primeseqs' );
	my $fiveargs = "-qspeedup 2 -dbcluster 80 0 $tir5name -p -d -seedlength 30 ";
	$fiveargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 10000 -identity 80";
	$vmatch_args{fivetir} = { seqs => \@fivetirs, args => $fiveargs };

	my $tir3name  = File::Spec->catfile( abs_path($dir), 'dbcluster-3primeseqs' );
	my $threeargs = "-qspeedup 2 -dbcluster 80 0 $tir3name -p -d -seedlength 30 ";
	$threeargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 10000 -identity 80";
	$vmatch_args{threetir} = { seqs => \@threetirs, args => $threeargs };

	# pdoms
	my $pdomname = File::Spec->catfile( abs_path($dir), 'dbcluster-pdoms' );
	my $pdomargs = "-qspeedup 2 -dbcluster 80 0 $pdomname -p -d -seedlength 30 -exdrop 3 ";
	$pdomargs .= "-l 40 -showdesc 0 -sort ld -best 10000";
	$vmatch_args{pdoms} = { seqs => \@pdoms, args => $pdomargs };
    }

    return \%vmatch_args;
}

sub cluster_features {
    my $self = shift;
    my ($dir) = @_;
    my $threads = $self->threads;

    my $args = $self->collect_feature_args($dir);
    $self->_remove_singletons($args);

    my $t0 = gettimeofday();
    my $doms = 0;
    my %reports;
    my $outfile = File::Spec->catfile( abs_path($dir), 'all_vmatch_reports.txt' );
    my $logfile = File::Spec->catfile( abs_path($dir), 'all_vmatch_reports.log' );
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $log, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";
    
    my $thr; 
    if ($threads % 3 == 0) {
	$thr = sprintf("%.0f", $threads/3);
    }
    elsif (+($threads-1) % 3 == 0) {
	$thr = sprintf("%.0f", $threads-1/3);
    }
    else {
	$thr = 1;
    }

    my $pm = Parallel::ForkManager->new($thr);
    local $SIG{INT} = sub {
        $log->warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
                              for my $bl (sort keys %$data_ref) {
                                  open my $report, '<', $bl or die "\n[ERROR]: Could not open file: $bl\n";
                                  print $out $_ while <$report>;
                                  close $report;
                                  unlink $bl;
                              }
                              my $t1 = gettimeofday();
                              my $elapsed = $t1 - $t0;
                              my $time = sprintf("%.2f",$elapsed/60);
                              say $log basename($ident),
                              " just finished with PID $pid and exit code: $exit_code in $time minutes";
                        } );

    for my $type (keys %$args) {
        for my $db (@{$args->{$type}{seqs}}) {
            $doms++;
            $pm->start($db) and next;
            $SIG{INT} = sub { $pm->finish };
            my $vmrep = $self->process_cluster_args($args, $type, $db);
            $reports{$vmrep} = 1;

            $pm->finish(0, \%reports);
        }
    }

    $pm->wait_all_children;
    close $out;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished running vmatch on $doms domains in $final_time minutes";
    close $log;

    return $outfile;
}

sub process_cluster_args {
    my $self = shift;
    my ($args, $type, $db) = @_;

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $index = File::Spec->catfile( abs_path($path), $name.'.index' );
    my $vmrep = File::Spec->catfile( abs_path($path), $name.'_vmatch-out.txt' );
    my $log   = File::Spec->catfile( abs_path($path), $name.'_vmatch-out.log' );

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($vmatchbin) = @{$config}{qw(vmatchbin)};
    my $vmatch  = File::Spec->catfile($vmatchbin, 'vmatch');
    my $mkvtree = File::Spec->catfile($vmatchbin, 'mkvtree');

    my $mkvtreecmd = "$mkvtree -db $db -dna -indexname $index -allout -v -pl ";
    if (defined $args->{$type}{prefixlen}) {
        $mkvtreecmd .= "$args->{$type}{prefixlen} ";
    }
    $mkvtreecmd .= "2>&1 > $log";
    say STDERR "DEBUG: $mkvtreecmd" if $self->debug;
    my $vmatchcmd  = "$vmatch $args->{$type}{args} $index > $vmrep";
    say STDERR "DEBUG: $vmatchcmd" if $self->debug;
    $self->run_cmd($mkvtreecmd);
    $self->run_cmd($vmatchcmd);
    unlink glob "$index*";
    unlink glob "$path/*.match";
    unlink $log;

    return $vmrep;
}

sub parse_clusters {
    my $self = shift;
    my ($clsfile) = @_;
    my $genome = $self->genome->absolute->resolve;

    # If there are no clusters, return here instead of processing the file.
    # Also, clean up the empty cluster and log file.
    unless (-s $clsfile) {
	my $logfile = $clsfile =~ s/.txt/.log/r;
	unlink $clsfile, $logfile;
	return {};
    }

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
        $name =~ s/$1//;
    }

    my (%cls, %all_seqs, %all_pdoms, $clusnum, $dom);
    open my $in, '<', $clsfile or die "\n[ERROR]: Could not open file: $clsfile\n";

    while (my $line = <$in>) {
        chomp $line;
        if ($line =~ /^# args=/) {
            my ($type) = ($line =~ /\/(\S+).index\z/);
            $dom = basename($type);
            $dom =~ s/${name}_//;
	    $dom =~ s/_pdom//;
	    $all_pdoms{$dom} = 1;
        }
        next if $line =~ /^# \d+/;
	if ($line =~ /^(\d+):/) {
	    $clusnum = $1;
        }
        elsif ($line =~ /^\s+(\S+)/) {
            my $element = $1;
            #$element =~ s/_\d+-?_?\d+$//;
	    # Below is a protein domain that was potentially concatenated from
	    # multiple domains of the same type. The coordinates of each
	    # subdomain are removed so the element and chromosome are left.
	    $element =~ s/_\d+-?_?\d+$//
		while $element =~ /_\d+-?_?\d+$/;
            push @{$cls{$dom}{$clusnum}}, $element;
        }
    }
    close $in;

    my (%elem_sorted, %multi_cluster_elem);
    for my $pdom (keys %cls) {
        for my $clsnum (keys %{$cls{$pdom}}) {
            for my $elem (@{$cls{$pdom}{$clsnum}}) {
                push @{$elem_sorted{$elem}}, { $pdom => $clsnum };
            }
        }
    }

    my %dom_orgs;
    for my $element (keys %elem_sorted) {
        my $string;
        my %pdomh;
        @pdomh{keys %$_} = values %$_ for @{$elem_sorted{$element}};
        for my $pdom (nsort keys %cls) {
            if (exists $pdomh{$pdom}) {
                $string .= length($string) ? "|$pdomh{$pdom}" : $pdomh{$pdom};
            }
            else {
                $string .= length($string) ? "|N" : "N";
            }
        }
        push @{$dom_orgs{$string}}, $element;
        undef $string;
    }

    #dd \%dom_orgs and exit;
    return \%dom_orgs;
}

sub _remove_singletons {
    my $self = shift;
    my ($args) = @_;

    my @singles;
    my ($index, $seqct) = (0, 0);
    for my $type (keys %$args) {
        delete $args->{$type} if ! @{$args->{$type}{seqs}};
        for my $db (@{$args->{$type}{seqs}}) {
            my $kseq = Bio::DB::HTS::Kseq->new($db);
            my $iter = $kseq->iterator();
            while (my $seqobj = $iter->next_seq) { $seqct++ if defined $seqobj->seq; }
            if ($seqct < 2) {
                push @singles, $index;
                #unlink $db; 
                # Do not remove the file of singletons. This is so we can keep track of all elements
		# for the solo-LTR search.
            }
            $index++;
            $seqct = 0;
        }

        if (@{$args->{$type}{seqs}}) {
            if (@singles > 1) {
                for (@singles) {
                    splice @{$args->{$type}{seqs}}, $_, 1;
                    @singles = map { $_ - 1 } @singles; # array length is changing after splice so we need to adjust offsets
                }
            }
            else {
                splice @{$args->{$type}{seqs}}, $_, 1 for @singles;
            }
        }
        else {
            delete $args->{$type};
        }
        $index = 0;
        @singles = ();
    }

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

    perldoc Tephra::Classify::Fams::Cluster


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
