package Tephra::Classify::Fams;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Sort::Naturally;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::DB::HTS::Kseq;
use Bio::DB::HTS::Faidx;
use Bio::GFF3::LowLevel qw(gff3_format_feature);
use List::Util          qw(min max);
use List::UtilsBy       qw(nsort_by);
use Time::HiRes         qw(gettimeofday);
use File::Path          qw(make_path);
use Cwd                 qw(abs_path);         
use Parallel::ForkManager;
use Carp 'croak';
use Try::Tiny;
use Tephra::Annotation::MakeExemplars;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Util',
     'Tephra::Role::GFF',
     'Tephra::Classify::Role::LogResults',
     'Tephra::Classify::Fams::Cluster',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Classify::Fams - Classify LTR/TIR transposons into families

=head1 VERSION

Version 0.10.0

=cut

our $VERSION = '0.10.0';
$VERSION = eval $VERSION;

=head1 SYNOPSIS

    use Tephra::Classify::Fams;

    my $genome  = 'genome.fasta';     # genome sequences in FASTA format
    my $outdir  = 'ltr_families_out'; # directory to place the results
    my $threads = 12;                 # the number of threads to use for parallel processing

    my $classify_fams_obj = Tephra::Classify::Fams->new(
        genome   => $genome,
        outdir   => $outdir,
        threads  => $threads,
    );

=cut

has genome => (
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

has type => (
      is       => 'ro',
      isa      => 'Str',
      required => 1,
      default  => 'LTR',
);

#
# methods
#
sub make_families {
    my $self = shift;
    my ($gff_obj, $log) = @_;

    my $outdir  = $self->outdir->absolute->resolve;
    my $threads = $self->threads;    
    my $tetype  = $self->type;
    #my $log     = $self->get_logger($logfile);
    
    my $t0 = gettimeofday();
    my $logfile = File::Spec->catfile($outdir, $tetype.'_superfamilies_thread_report.log');
    open my $fmlog, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";

    my $pm = Parallel::ForkManager->new(3);
    local $SIG{INT} = sub {
        $log->warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    my (%reports, @family_fastas, @annotated_ids);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      for my $type (keys %$data_ref) {
				  my $family_stats = $data_ref->{$type}{family_stats};
				  my ($sf, $elemct, $famct, $famtot, $singct) =
                                      @{$family_stats}{qw(superfamily total_elements families total_in_families singletons)};
				  my ($sfam) = ($sf =~ /_?((?:\w+\d+\-)?\w+)\z/);
				  $sfam =~ s/.*hat.*/hAT/i;
				  $sfam =~ s/.*tc1-mariner.*/Tc1-Mariner/i;
				  $sfam =~ s/.*mutator.*/Mutator/i;
				  $sfam =~ s/.*cacta.*/CACTA/i;
				  $sfam =~ s/.*gypsy.*/Gypsy/i;
				  $sfam =~ s/.*copia.*/Copia/i;
				  $sfam =~ s/.*unclassified.*/unclassified/i;
				  $sfam =~ s/.*mite.*/MITE/i;

				  my $tot_str = sprintf("%-70s %-10s", "Results - Number of $sfam families:", $famct);
				  my $fam_str = sprintf("%-70s %-10s", "Results - Number of $sfam elements in families:", $famtot);
				  my $sng_str = sprintf("%-70s %-10s", "Results - Number of $sfam singleton families/elements:", $singct);
				  my $num_str = sprintf("%-70s %-10s", "Results - Number of $sfam elements (for debugging):", $elemct);

				  $log->info($tot_str);
				  $log->info($fam_str);
				  $log->info($sng_str);
				  $log->info($num_str);

				  push @family_fastas, $data_ref->{$type}{family_fasta};
				  push @annotated_ids, $data_ref->{$type}{annotated_ids};
			      }
			      my $t1 = gettimeofday();
			      my $elapsed = $t1 - $t0;
			      my $time = sprintf("%.2f",$elapsed/60);
			      say $fmlog "$ident just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    for my $type (keys %$gff_obj) {
	$pm->start($type) and next;
	$SIG{INT} = sub { $pm->finish };

	my ($fams, $ids, $family_stats) = 
	    $self->run_family_classification($gff_obj->{$type}, $tetype);

	$reports{$type} = { family_fasta   => $fams, 
			    annotated_ids  => $ids, 
	                    family_stats   => $family_stats };

	$pm->finish(0, \%reports);
    }

    $pm->wait_all_children;

    my (%outfiles, %annot_ids);
    @outfiles{keys %$_}  = values %$_ for @family_fastas;
    @annot_ids{keys %$_} = values %$_ for @annotated_ids;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $fmlog "\n========> Finished classifying $tetype families in $final_time minutes";
    $log->info("Finished classifying $tetype families in $final_time minutes.");
    close $fmlog;

    return (\%outfiles, \%annot_ids);
}

sub run_family_classification {
    my $self = shift;
    my ($gff, $tetype) = @_;
    my $genome  = $self->genome->absolute->resolve;
    my $threads = $self->threads;
    my $debug   = $self->debug;

    my ($dir, $pdom_famid_map);
    if ($tetype eq 'LTR') {
	($dir, $pdom_famid_map) = $self->extract_ltr_features($gff);
    }
    if ($tetype eq 'TIR') {
	($dir, $pdom_famid_map) = $self->extract_tir_features($gff);
    }
    #dd $pdom_famid_map and exit;

    my $clusters = $self->cluster_features($dir);
    my $dom_orgs = $self->parse_clusters($clusters);
    my ($dom_fam_map, $fas_obj, $unmerged_stats) = 
	$self->make_fasta_from_dom_orgs($dom_orgs, $clusters, $tetype);
    
    my $blastout = $self->process_blast_args($fas_obj);
    my $matches  = $self->parse_blast($blastout);

    my ($fams, $ids, $merged_stats) = 
	$self->write_families($pdom_famid_map, $dom_fam_map, $unmerged_stats, $matches, $clusters, $tetype);

    my $exm_obj = Tephra::Annotation::MakeExemplars->new(
        genome  => $genome,
        dir     => $dir,
        gff     => $gff,
	threads => $threads,
	debug   => $debug,
    );

    $exm_obj->make_exemplars;

    my %families = (
	$fas_obj->{family_fasta}    => 1,
	$fas_obj->{singleton_fasta} => 1,
    );

    return (\%families, $ids, $merged_stats);
}

sub make_fasta_from_dom_orgs {
    my $self = shift;
    my ($dom_orgs, $clsfile, $tetype) = @_;

    my ($cname, $cpath, $csuffix) = fileparse($clsfile, qr/\.[^.]*/);
    my $dir  = basename($cpath);
    my ($sf) = ($dir =~ /_((?:\w+\d+\-)?\w+)$/);
    unless (defined $sf) {
	say STDERR "\n[ERROR]: Can not get superfamily from $clsfile at $.";
    }

    my $sfname;
    if ($tetype eq 'LTR') {
	$sfname = 'RLG' if $sf =~ /gypsy/i;
	$sfname = 'RLC' if $sf =~ /copia/i;
	$sfname = 'RLX' if $sf =~ /unclassified/i;
    }
    if ($tetype eq 'TIR') {
	$sfname = 'DTA' if $sf =~ /hat/i;
	$sfname = 'DTC' if $sf =~ /cacta/i;
	$sfname = 'DTM' if $sf =~ /mutator/i;
	$sfname = 'DTT' if $sf =~ /mariner/i;
	$sfname = 'DTX' if $sf =~ /unclassified/i;
	$sfname = 'DTX' if $sf =~ /mite/i;
    }

    my @compfiles;
    find( sub { push @compfiles, $File::Find::name if /complete.fasta$/ }, $cpath );
    my $ltrfas = shift @compfiles;
    my $seqstore = $self->_store_seq($ltrfas);
    my $elemct = (keys %$seqstore);

    my $famfile = $sfname.'_families.fasta';
    my $foutfile = File::Spec->catfile($cpath, $famfile);
    open my $out, '>>', $foutfile or die "\n[ERROR]: Could not open file: $foutfile\n";

    my %fam_map;
    my ($idx, $famtot) = (0, 0);
    for my $str (reverse sort { @{$dom_orgs->{$a}} <=> @{$dom_orgs->{$b}} } keys %$dom_orgs) {  
	my $famnum = $sfname."_family$idx";
        for my $elem (@{$dom_orgs->{$str}}) {
            if (exists $seqstore->{$elem}) {
                $famtot++;
                my $coordsh = $seqstore->{$elem};
                my $coords  = (keys %$coordsh)[0];
                $seqstore->{$elem}{$coords} =~ s/.{60}\K/\n/g;
                say $out join "\n", ">$sfname"."_family$idx"."_$elem"."_$coords", $seqstore->{$elem}{$coords};
                delete $seqstore->{$elem};
		push @{$fam_map{$famnum}}, $elem;
            }
            else {
                croak "\n[ERROR]: $elem not found in store. Exiting.";
            }
        }
        $idx++;
    }
    close $out;
    my $famct = $idx;
    $idx = 0;

    my $reduc = (keys %$seqstore);
    my $singfile = $sfname.'_singletons.fasta';
    my $soutfile = File::Spec->catfile( abs_path($cpath), $singfile );
    open my $outx, '>>', $soutfile or die "\n[ERROR]: Could not open file: $soutfile\n";

    if (%$seqstore) {
        for my $k (nsort keys %$seqstore) {
            my $coordsh = $seqstore->{$k};
            my $coords  = (keys %$coordsh)[0];
            $seqstore->{$k}{$coords} =~ s/.{60}\K/\n/g;
            say $outx join "\n", ">$sfname"."_singleton_family$idx"."_$k"."_$coords", $seqstore->{$k}{$coords};
            $idx++;
        }
    }
    close $outx;
    my $singct = $idx;
    undef $seqstore;

    return (\%fam_map,
	    { family_fasta      => $foutfile, 
	      singleton_fasta   => $soutfile, 
	      family_count      => $famct, 
	      total_in_families => $famtot, 
	      singleton_count   => $singct }, 
	   { superfamily       => $sf,
	     total_elements    => $elemct,
	     families          => $famct,
	     total_in_families => $famtot,
	     singletons        => $singct });
}

sub process_blast_args {
    my $self = shift;
    my ($obj) = @_;
    my $threads = $self->threads;
    my ($query, $db) = @{$obj}{qw(singleton_fasta family_fasta)};
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

    my $blast_report = $self->run_blast({ query => $query, db => $blastdb, threads => $thr, outfile => $report, sort => 'bitscore' });
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
    open my $in, '<', $blast_report or die "\n[ERROR]: Could not open file: $blast_report\n";
    while (my $line = <$in>) {
	chomp $line;
	my ($queryid, $hitid, $pid, $hitlen, $mmatchct, $gapct, 
	    $qhit_start, $qhit_end, $hhit_start, $hhit_end, $evalue, $score) = split /\t/, $line;
	my ($qstart, $qend) = ($queryid =~ /(\d+)-?_?(\d+)$/);
	my ($hstart, $hend) = ($hitid =~ /(\d+)-?_?(\d+)$/);
	my $qlen = $qend - $qstart + 1;
	my $hlen = $hend - $hstart + 1;
	my $minlen = min($qlen, $hlen); # we want to measure the coverage of the smaller element
	my ($coords) = ($queryid =~ /_(\d+_\d+)$/);
        $queryid =~ s/_$coords//;
	my ($family) = ($hitid =~ /(\w{3}_(?:singleton_)?family\d+)_/);
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
    my ($pdom_famid_map, $dom_fam_map, $unmerged_stats, $matches, $clsfile, $tetype) = @_;
    my $gff = $self->gff;

    #dd $pdom_famid_map and exit;
    my ($cname, $cpath, $csuffix) = fileparse($clsfile, qr/\.[^.]*/);
    my $dir = basename($cpath);
    #my ($sf) = ($dir =~ /_(\w+)$/);
    my ($sf) = ($dir =~ /_((?:\w+\d+\-)?\w+)$/);
    unless (defined $sf) {
        say STDERR "\n[ERROR]: Can not get superfamily name from $clsfile at $.";
    }

    my $sfname;
    if ($tetype eq 'LTR') {
	$sfname = 'RLG' if $sf =~ /gypsy/i;
	$sfname = 'RLC' if $sf =~ /copia/i;
	$sfname = 'RLX' if $sf =~ /unclassified/i;
    }
    if ($tetype eq 'TIR') {
	$sfname = 'DTA' if $sf =~ /hat/i;
	$sfname = 'DTC' if $sf =~ /cacta/i;
	$sfname = 'DTM' if $sf =~ /mutator/i;
	$sfname = 'DTT' if $sf =~ /mariner/i;
	$sfname = 'DTX' if $sf =~ /unclassified/i;
	$sfname = 'DTX' if $sf =~ /mite/i;
    }

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $domoutfile = File::Spec->catfile($path, $name.'_family-level_domain_org.tsv');
    open my $domf, '>>', $domoutfile or die "\n[ERROR]: Could not open file: $domoutfile\n";

    my @compfiles;
    find( sub { push @compfiles, $File::Find::name if /complete.fasta$/ }, $cpath );
    my $ltrfas = shift @compfiles;
    my $seqstore = $self->_store_seq($ltrfas);
    my $elemct = (keys %$seqstore);
	
    my ($idx, $famtot, $tomerge) = (0, 0, 0);
    my (%fastas, %annot_ids);

    my $fam_id_map;
    if (defined $matches) {
	my $tomerge = $self->_compare_merged_nonmerged($matches, $seqstore, $unmerged_stats);
	$fam_id_map = $tomerge == 1 ? $matches : $dom_fam_map;
    }

    my $re = qr/helitron\d+|non_LTR_retrotransposon\d+|(?:LTR|TRIM|LARD)_retrotransposon\d+|terminal_inverted_repeat_element\d+|MITE\d+/;
    #my ($element) = ($elemnum =~ /($re)/);
    #say "ltrfas: $ltrfas";
    #say join "\n", keys %$seqstore and exit;
    #dd $fam_id_map and exit;
    #for my $str (reverse sort { @{$fam_id_map->{$a}} <=> @{$fam_id_map->{$b}} } keys %$fam_id_map) {
    for my $str (reverse nsort keys %$fam_id_map) {
	my $famfile = $sf."_family$idx".".fasta";
	my $outfile = File::Spec->catfile( abs_path($cpath), $famfile );
	open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	for my $elem (@{$fam_id_map->{$str}}) {
	    my $query = $elem =~ s/\w{3}_singleton_family\d+_//r;
	    if (exists $seqstore->{$query}) {
		$famtot++;
		my $coordsh = $seqstore->{$query};
		my $coords  = (keys %$coordsh)[0];
		$seqstore->{$query}{$coords} =~ s/.{60}\K/\n/g;
		say $out join "\n", ">$sfname"."_family$idx"."_$query"."_$coords", $seqstore->{$query}{$coords};
		delete $seqstore->{$query};
		$annot_ids{$query} = $sfname."_family$idx";
		my ($element) = ($query =~ /($re)/);

		if (exists $pdom_famid_map->{$element}) { 
		    $self->write_fam_pdom_organization({ famid_map => $pdom_famid_map, 
							 outfh     => $domf, 
							 elemnum   => $query, 
							 elemid    => "$sfname"."_family$idx"."_$query"."_$coords", 
							 famname   => $sfname."_family$idx"});
		}
	    }
	    else {
		croak "\n[ERROR]: $query not found in store. Exiting.";
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
	my $xoutfile = File::Spec->catfile( abs_path($cpath), $famxfile );
	open my $outx, '>', $xoutfile or die "\n[ERROR]: Could not open file: $xoutfile\n";
	for my $k (nsort keys %$seqstore) {
	    my $coordsh = $seqstore->{$k};
	    my $coords  = (keys %$coordsh)[0];
	    $seqstore->{$k}{$coords} =~ s/.{60}\K/\n/g;
	    say $outx join "\n", ">$sfname"."_singleton_family$idx"."_$k"."_$coords", $seqstore->{$k}{$coords};
	    $annot_ids{$k} = $sfname."_singleton_family$idx";
	    my ($element) = ($k =~ /($re)/);

	    if (exists $pdom_famid_map->{$element}) {
		$self->write_fam_pdom_organization({ famid_map => $pdom_famid_map, 
						     outfh     => $domf, 
						     elemnum   => $k,
						     elemid    => "$sfname"."_singleton_family$idx"."_$k"."_$coords", 
						     famname   => $sfname."_singleton_family$idx" });
	    }
	    $idx++;
	}
	close $outx;
	$fastas{$xoutfile} = 1;
    }
    my $singct = $idx;
    close $domf;
    unlink $domoutfile unless -s $domoutfile;

    return (\%fastas, \%annot_ids,
	    { superfamily       => $sf,
	      total_elements    => $elemct,
	      families          => $famct,
	      total_in_families => $famtot,
	      singletons        => $singct });
}

sub combine_families {
    my ($self) = shift;
    my ($outfiles) = @_;
    my $genome = $self->genome->absolute->resolve;
    my $outdir = $self->outdir->absolute->resolve;
    my $outgff = $self->gff;
    
    my ($name, $path, $suffix) = fileparse($outgff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'.fasta');

    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";

    for my $file (nsort keys %$outfiles) {
	unlink $file && next unless -s $file;
	my $kseq = Bio::DB::HTS::Kseq->new($file);
	my $iter = $kseq->iterator();
	while (my $seqobj = $iter->next_seq) {
	    my $id  = $seqobj->name;
	    my $seq = $seqobj->seq;
	    $seq =~ s/.{60}\K/\n/g;
	    say $out join "\n", ">$id", $seq;
	}
	#unlink $file;
    }
    close $outfile;

    return;
}

sub annotate_gff {
    my $self = shift;
    my ($annot_ids, $index, $ingff, $type) = @_;
    my $outdir = $self->outdir->absolute->resolve;
    my $outgff = $self->gff;

    my $new_type;
    if ($type eq 'TIR') {
	$new_type = 'MITE';

    }
    elsif ($type eq 'LTR') { 
	$new_type = 'LARD_retrotransposon';
    }
    else {
	say "\n[ERROR]: Could not determine LTR/TIR type to annotate GFF3. This is a bug, please report it. Exiting.\n";
	exit(1);

    }

    my ($header, $features) = $self->collect_gff_features($ingff);
    #open my $in, '<', $ingff or die "\n[ERROR]: Could not open file: $ingff\n";
    open my $out, '>', $outgff or die "\n[ERROR]: Could not open file: $outgff\n";
    say $out $header;

    my $is_lard_mite = 0;
    my ($new_id, $gff_str, $seq_id, $strand, $source);
    #dd $features and exit;
    for my $rep_region (nsort_by { m/repeat_region\d+\|\|(\d+)\|\|\d+/ and $1 } keys %$features) {
	my ($rreg_id, $s, $e) = split /\|\|/, $rep_region;
        for my $feature (@{$features->{$rep_region}}) {
	    if ($feature->{type} =~ /(?:LTR|TRIM)_retrotransposon|terminal_inverted_repeat_element|MITE/) {
                #my ($seq_id, $start, $end) = @{$ltr_feature}{qw(seq_id start end)};
		#dd $feature and exit;
		$is_lard_mite = 1;
		
                my $id = $feature->{attributes}{ID}[0];
		($seq_id, $strand, $source) = @{$feature}{qw(seq_id strand source)};

		if (exists $index->{$id}) {
		    $feature->{type} = $new_type;
                    $new_id = $index->{$id};
		}
		else {
		    $new_id = $id;
		}

		my $key = join "_", $new_id, $seq_id;
		if (exists $annot_ids->{$key}) {
		    $feature->{attributes}{ID}[0] = $new_id;
		    $feature->{attributes}{family}[0] = $annot_ids->{$key};		    
		}
	    }
	    else {
		if ($is_lard_mite) {
		    $feature->{attributes}{Parent}[0] = $new_id;
		}
	    }
	    my $gff_feat = gff3_format_feature($feature);
	    $gff_str .= $gff_feat;
	    #chomp $gff_str;
	    #say $out $gff_str;
	}
	chomp $gff_str;
	say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg_id";
	say $out $gff_str;
	$is_lard_mite = 0;
	undef $gff_str;
    }
    ##debug
	    #my $is_lard = 0;
    #my $new_id;
    #dd $annot_ids;
    #while (my $line = <$in>) {
	#chomp $line;
	#if ($line =~ /^#/) {
	    #say $out $line;
	#}
	#else {
	    #my @f = split /\t/, $line;
	    #if ($f[8] =~ /((?:LTR|TRIM|LARD)_retrotransposon\d+|terminal_inverted_repeat_element\d+)/) {
		#$my ($id) = ($f[8] =~ /ID=((?:LTR|TRIM|LARD)_retrotransposon\d+|terminal_inverted_repeat_element\d+);/);
		#my $key  = $id."_$f[0]";
		#my $id = $1;
		
		#if (exists $lard_index->{$id}) {
		    #$new_id = $lard_index->{$id};
		    #say "found in lard index: $id";
		    #$is_lard = 1;
		    #$f[2] =~ s/LTR/LARD/;
		    #$f[8] =~ s/LTR/LARD/;
		#}
		#else {
		#    $new_id = $id;
		#}

		#my $key = $new_id."_$f[0]";
		#say join q{ }, $key, $id if $is_lard;
		#if (exists $annot_ids->{$key}) {
		    #my $family = $annot_ids->{$key};
		    #$f[8] =~ s/ID=$id\;/ID=$new_id;family=$family;/;
		    #$f[8] =~ s/Parent=$id\;/Parent=$new_id/;
		    #say $out join "\t", @f;
		#}
		#else {
		#    say $out join "\t", @f;
		#}
	    #}
	    #else {
		#say $out join "\t", @f;
	    #}
	    #$is_lard = 0;
	#}
	
    #}
    #close $in;
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
		croak "\n[ERROR]: $query not found in store. Exiting.";
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

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Classify::Fams


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 
	
__PACKAGE__->meta->make_immutable;

1;
