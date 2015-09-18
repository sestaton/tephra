package Tephra::Classify::LTRFams;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use Time::HiRes qw(gettimeofday);
use Parallel::ForkManager;
use Cwd;
use Try::Tiny;
use namespace::autoclean;

use Data::Dump;
use Data::Printer;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

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

#
# methods
#
sub extract_features {
    my $self = shift;
    my $fasta  = $self->genome;
    my $dir    = $self->outdir;
    my ($infile) = @_;
    
    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $comp = File::Spec->catfile($dir, $name."_complete.fasta");
    my $ppts = File::Spec->catfile($dir, $name."_ppt.fasta");
    my $pbs  = File::Spec->catfile($dir, $name."_pbs.fasta");
    my $five_pr_ltrs  = File::Spec->catfile($dir, $name."_5prime-ltrs.fasta");
    my $three_pr_ltrs = File::Spec->catfile($dir, $name."_3prime-ltrs.fasta");

    open my $allfh, '>>', $comp;
    open my $pptfh, '>>', $ppts;
    open my $pbsfh, '>>', $pbs;
    open my $fivefh, '>>', $five_pr_ltrs;
    open my $threfh, '>>', $three_pr_ltrs;

    my $gffio = Bio::Tools::GFF->new( -file => $infile, -gff_version => 3 );

    my ($start, $end, $elem_id, $key, %feature, %ltrs, %seen);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'LTR_retrotransposon') {
	    my @string = split /\t/, $feature->gff_string;
	    ($elem_id) = ($string[8] =~ /ID=?\s+?(LTR_retrotransposon\d+)/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $key = join ".", $elem_id, $start, $end;
	    $ltrs{$key}{'full'} = join "-", $string[0], $feature->primary_tag, @string[3..4];
	}
	next unless defined $start && defined $end;
	if ($feature->primary_tag eq 'long_terminal_repeat') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my $ltrkey = join "||", $string[0], $feature->primary_tag, @string[3..4];
		push @{$ltrs{$key}{'ltrs'}}, $ltrkey unless exists $seen{$ltrkey};
		$seen{$ltrkey} = 1;
	    }
	}
	elsif ($feature->primary_tag eq 'primer_binding_site') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my ($name) = ($string[8] =~ /trna \"?(\w+.*)\"?\s+\;/);
		$ltrs{$key}{'pbs'} =
		    join "||", $string[0], $feature->primary_tag, $name, @string[3..4];
	    }
	}
	elsif ($feature->primary_tag eq 'protein_match') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my ($name) = ($string[8] =~ /name \"?(\w+)\"?/);
		my $pdomkey = join "||", $string[0], $feature->primary_tag, $name, @string[3..4];
		push @{$ltrs{$key}{'pdoms'}}, $pdomkey unless exists $seen{$pdomkey};
		$seen{$pdomkey} = 1;
	    }
	}
	elsif ($feature->primary_tag eq 'RR_tract') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		$ltrs{$key}{'ppt'} =
		    join "||", $string[0], $feature->primary_tag, @string[3..4];
	    }
	}
    }

    #dd \%ltrs;
    
    my %pdoms;
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
	my ($element, $rstart, $rend) = split /\./, $ltr;
	# full element
	my ($source, $prim_tag, $start, $end) = split /\-/, $ltrs{$ltr}{'full'};
	my $outfile = File::Spec->catfile($dir, $ltr.".fasta");
	$self->subseq($fasta, $source, $element, $start, $end, $outfile, $allfh);

	# pbs
	if ($ltrs{$ltr}{'pbs'}) {
	    my ($pbssource, $pbstag, $trna, $pbsstart, $pbsend) = split /\|\|/, $ltrs{$ltr}{'pbs'};
	    my $pbs_tmp = File::Spec->catfile($dir, $ltr."_pbs.fasta");
	    $self->subseq($fasta, $pbssource, $element, $pbsstart, $pbsend, $pbs_tmp, $pbsfh);
	}

	# ppt
	if ($ltrs{$ltr}{'ppt'}) {
	    my ($pptsource, $ppttag, $pptstart, $pptend) = split /\|\|/, $ltrs{$ltr}{'ppt'};
	    my $ppt_tmp = File::Spec->catfile($dir, $ltr."_ppt.fasta");
	    $self->subseq($fasta, $source, $element, $pptstart, $pptend, $ppt_tmp, $pptfh);
	}

	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltrtag, $s, $e) = split /\|\|/, $ltr_repeat;
	    if ($ltrct) {
		my $fiveprime_tmp = File::Spec->catfile($dir, $ltr."_5prime-ltr.fasta");
		$self->subseq($fasta, $src, $element, $s, $e, $fiveprime_tmp, $fivefh);
	    }
	    else {
		my $threeprime_tmp = File::Spec->catfile($dir, $ltr."_3prime-ltr.fasta");
		$self->subseq($fasta, $src, $element, $s, $e, $threeprime_tmp, $threfh);
		$ltrct++;
	    }
	}
	$ltrct = 0;

	if ($ltrs{$ltr}{'pdoms'}) {
	    for my $ltr_repeat (@{$ltrs{$ltr}{'pdoms'}}) {
		my ($src, $pdomtag, $name, $s, $e ) = split /\|\|/, $ltr_repeat;
		#"Ha10||protein_match||UBN2||132013916||132014240",
		push @{$pdoms{$name}}, join "||", $src, $element, $s, $e;
	    }
	}
    }
    close $allfh;
    close $pptfh;
    close $pbsfh;
    close $fivefh;
    close $threfh;

    for my $pdom_type (keys %pdoms) {
	my $pdom_file = File::Spec->catfile($dir, $pdom_type."_pdom.fasta");
	open my $fh, '>>', $pdom_file;
	for my $ltrpdom (@{$pdoms{$pdom_type}}) {
	    my ($src, $elem, $s, $e) = split /\|\|/, $ltrpdom;
	    my $tmp = File::Spec->catfile($dir, $elem."_".$pdom_type.".fasta");
	    $self->subseq($fasta, $src, $elem, $s, $e, $tmp, $fh);
	}
	close $fh;
	unlink $pdom_file if ! -s $pdom_file;
    }

    for my $file ($comp, $ppts, $pbs, $five_pr_ltrs, $three_pr_ltrs) {
	unlink $file if ! -s $file;
    }
}

sub collect_feature_args {
    my $self = shift;
    my $dir = $self->outdir;
    my (@fiveltrs, @threeltrs, @ppt, @pbs, @pdoms, %vmatch_args);
    find( sub { push @fiveltrs, $File::Find::name if -f and /5prime-ltrs.fasta$/ }, $dir);
    find( sub { push @threeltrs, $File::Find::name if -f and /3prime-ltrs.fasta$/ }, $dir);
    find( sub { push @ppt, $File::Find::name if -f and /ppts.fasta$/ }, $dir);
    find( sub { push @pbs, $File::Find::name if -f and /pbs.fasta$/ }, $dir);
    find( sub { push @pdoms, $File::Find::name if -f and /pdom.fasta$/ }, $dir);

    # ltr
    my $ltr5name = File::Spec->catfile($dir, 'dbcluster-5primeseqs');
    my $fiveargs = "-dbcluster 80 20 $ltr5name -p -d -seedlength 10 ";
    $fiveargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 100 -identity 80";
    $vmatch_args{fiveltr} = { seqs => \@fiveltrs, args => $fiveargs };

    my $ltr3name  = File::Spec->catfile($dir, 'dbcluster-3primeseqs');
    my $threeargs = "-dbcluster 80 20 $ltr3name -p -d -seedlength 10 ";
    $threeargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 100 -identity 80";
    $vmatch_args{threeltr} = { seqs => \@threeltrs, args => $threeargs };

    # pbs/ppt
    my $pbsname = File::Spec->catfile($dir, 'dbcluster-pbs');
    my $pbsargs = "-dbcluster 90 90 $pbsname -p -d -seedlength 5 -exdrop 2 ";
    $pbsargs .= "-l 3 -showdesc 0 -sort ld -best 100";
    $vmatch_args{pbs} = { seqs => \@pbs, args => $pbsargs, prefixlen => 3 };

    my $pptname = File::Spec->catfile($dir, 'dbcluster-ppt');
    my $pptargs = "-dbcluster 90 90 $pptname -p -d -seedlength 5 -exdrop 2 ";
    $pptargs .= "-l 3 -showdesc 0 -sort ld -best 100";
    $vmatch_args{ppt} = { seqs => \@ppt, args => $pptargs, prefixlen => 5 };

    # pdoms
    my $pdomname = File::Spec->catfile($dir, 'dbcluster-pdoms');
    my $pdomargs = "-dbcluster 80 80 $pdomname -p -d -seedlength 10 -exdrop 3 ";
    $pdomargs .= "-l 40 -showdesc 0 -sort ld -best 100";
    $vmatch_args{pdoms} = { seqs => \@pdoms, args => $pdomargs };

    return \%vmatch_args;
}

sub cluster_features {
    my $self = shift;
    my $dir  = $self->outdir;
    my $threads = $self->threads;
    #my ($args) = @_;

    my $args = $self->collect_feature_args;
    $self->_remove_singletons($args);

    #use Data::Dump;
    #say STDERR "DEBUG: ";
    #dd $args;

    my $t0 = gettimeofday();
    my $doms = 0;
    my %reports;
    my $outfile = File::Spec->catfile($dir, 'all_vmatch_reports.txt');
    my $logfile = File::Spec->catfile($dir, 'all_vmatch_reports.log');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $log, '>>', $logfile or die "\nERROR: Could not open file: $logfile\n";
    
    my $pm = Parallel::ForkManager->new($threads);
    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      for my $bl (sort keys %$data_ref) {
				  open my $report, '<', $bl or die "\nERROR: Could not open file: $bl\n";
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
    my $index = File::Spec->catfile($path, $name.".index");
    my $vmrep = File::Spec->catfile($path, $name."_vmatch-out.txt");
    my $log   = File::Spec->catfile($path, $name."_vmatch-out.log");;

    my $mkvtreecmd = "mkvtree -db $db -dna -indexname $index -allout -v -pl ";
    if (defined $args->{$type}{prefixlen}) {
	$mkvtreecmd .= "$args->{$type}{prefixlen} ";
    }
    $mkvtreecmd .= "2>&1 > $log";
    my $vmatchcmd  = "vmatch $args->{$type}{args} $index > $vmrep";
    $self->run_cmd($mkvtreecmd);
    $self->run_cmd($vmatchcmd);
    unlink glob "$index*";

    return $vmrep;
}

sub subseq {
    my $self = shift;
    my ($fasta, $loc, $elem, $start, $end, $tmp, $out) = @_;
    my $cmd = "samtools faidx $fasta $loc:$start-$end > $tmp";
    $self->run_cmd($cmd);

    my $id = join "_", $loc, $elem, "$start-$end";
    if (-s $tmp) {
	my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
	while (my $seqobj = $seqio->next_seq) {
	    my $seq = $seqobj->seq;
	    if ($seq) {
		$seq =~ s/.{60}\K/\n/g;
		say $out join "\n", ">".$id, $seq;
	    }
	}
    }
    unlink $tmp;
}

sub _remove_singletons {
    my $self = shift;
    my ($args) = @_;

    my @singles;
    my ($index, $seqct) = (0, 0);
    for my $type (keys %$args) {
	delete $args->{$type} if ! @{$args->{$type}{seqs}};
	for my $db (@{$args->{$type}{seqs}}) {
	    my $seqio = Bio::SeqIO->new( -file => $db, -format => 'fasta' );
	    while (my $seqobj = $seqio->next_seq) { $seqct++ if defined $seqobj->seq; }
	    if ($seqct < 2) {
		push @singles, $index;
		#say STDERR "under 2: $db";
		unlink $db;
	    }
	    $index++;
	    $seqct = 0;
	}

	if (@{$args->{$type}{seqs}}) {
	    #p @{$args->{$type}{seqs}};
	    #p @singles;
	    if (@singles > 1) {
		for (@singles) {
		    splice @{$args->{$type}{seqs}}, $_, 1;
		    @singles = map { $_ - 1 } @singles; # array length is changing after splice so we need to adjust offsets
		}
	    }
	    else {
		splice @{$args->{$type}{seqs}}, $_, 1 for @singles;
	    }
	    #p @{$args->{$type}{seqs}};
	}
	else {
	    delete $args->{$type};
	}
	$index = 0;
	@singles = ();
    }
    #p $args;
}
	
__PACKAGE__->meta->make_immutable;

1;
