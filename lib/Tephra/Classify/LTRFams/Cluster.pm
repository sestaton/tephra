package Tephra::Classify::LTRFams::Cluster;

use 5.010;
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
use Parallel::ForkManager;
use Carp 'croak';
use Try::Tiny;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

=head1 NAME

     Tephra::Classify::LTRFams::Cluster - Helper role with clustering routines for LTR classification

=head1 VERSION

Version 0.04.3

=cut

our $VERSION = '0.04.3';
$VERSION = eval $VERSION;

sub extract_features {
    my $self = shift;
    my $fasta  = $self->genome;
    my $dir    = $self->outdir;
    my ($infile) = @_;
    
    my $index = $self->index_ref($fasta);

    my ($name, $path, $suffix) = fileparse($infile, qr/\.[^.]*/);
    my $type = ($name =~ /(?:gypsy|copia|unclassified)$/i);
    die "\nERROR: Unexpected input. Should match /gypsy|copia|unclassified$/i. Exiting."
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

    open my $allfh, '>>', $comp or die "\nERROR: Could not open file: $comp\n";
    open my $pptfh, '>>', $ppts or die "\nERROR: Could not open file: $ppts\n";
    open my $pbsfh, '>>', $pbs or die "\nERROR: Could not open file: $pbs\n";
    open my $fivefh, '>>', $five_pr_ltrs or die "\nERROR: Could not open file: $five_pr_ltrs\n";
    open my $threfh, '>>', $three_pr_ltrs or die "\nERROR: Could not open file: $three_pr_ltrs\n";

    open my $gffio, '<', $infile or die "\nERROR: Could not open file: $infile\n";

    my (%feature, %ltrs, %coord_map, %seen);
    while (my $line = <$gffio>) {
        chomp $line;
        next if $line =~ /^#/;
        my $feature = gff3_parse_feature($line);

        if ($feature->{type} eq 'LTR_retrotransposon') {
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
            my $name = $feature->{attributes}{trna};
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
                push @{$ltrs{$pkey}{'pdoms'}{$name}}, $pdomkey unless exists $seen{$pdomkey};
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

    my (%pdoms, %seen_pdoms);
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
        my ($element, $rstart, $rend) = split /\|\|/, $ltr;
        # full element
        my ($source, $prim_tag, $fstart, $fend) = split /\|\|/, $ltrs{$ltr}{'full'};
        $self->subseq($index, $source, $element, $fstart, $fend, $allfh);

        # pbs
        if ($ltrs{$ltr}{'pbs'}) {
            my ($pbssource, $pbstag, $trna, $pbsstart, $pbsend) = split /\|\|/, $ltrs{$ltr}{'pbs'};
            $self->subseq($index, $pbssource, $element, $pbsstart, $pbsend, $pbsfh);
        }

        # ppt
        if ($ltrs{$ltr}{'ppt'}) {
            my ($pptsource, $ppttag, $pptstart, $pptend) = split /\|\|/, $ltrs{$ltr}{'ppt'};
            $self->subseq($index, $source, $element, $pptstart, $pptend, $pptfh);
        }

        for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
            my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
            if ($ltrct) {
                $self->subseq($index, $src, $element, $s, $e, $fivefh);
                $ltrct = 0;
            }
            else {
                $self->subseq($index, $src, $element, $s, $e, $threfh);
                $ltrct++;
            }
        }

	if ($ltrs{$ltr}{'pdoms'}) {
            for my $model_name (keys %{$ltrs{$ltr}{'pdoms'}}) {
                for my $ltr_repeat (@{$ltrs{$ltr}{'pdoms'}{$model_name}}) {
                    my ($src, $pdomtag, $name, $s, $e, $str) = split /\|\|/, $ltr_repeat;
                    #"Ha10||protein_match||UBN2||132013916||132014240|+",
                    next if $model_name =~ /transpos(?:ase)?|mule|(?:dbd|dde)?_tnp_(?:hat)?|duf4216/i; 
                    # The above is so we do not classify elements based domains derived from or belonging to DNA transposons
                    push @{$pdoms{$src}{$element}{$model_name}}, join "||", $s, $e, $str;
                }
            }
        }
    }
    close $allfh;
    close $pptfh;
    close $pbsfh;
    close $fivefh;
    close $threfh;

    ## This is where we merge overlapping hits in a chain and concatenate non-overlapping hits
    ## to create a single domain sequence for each element
    for my $src (keys %pdoms) {
        for my $element (keys %{$pdoms{$src}}) {
            my ($pdom_s, $pdom_e, $str);
            for my $pdom_type (keys %{$pdoms{$src}{$element}}) {
                my (%lrange, %seqs, $union);
                my $pdom_file = File::Spec->catfile($resdir, $pdom_type.'_pdom.fasta');
                open my $fh, '>>', $pdom_file or die "\nERROR: Could not open file: $pdom_file\n";
                for my $split_dom (@{$pdoms{$src}{$element}{$pdom_type}}) {
                    ($pdom_s, $pdom_e, $str) = split /\|\|/, $split_dom;
                    push @{$lrange{$src}{$element}{$pdom_type}}, "$pdom_s..$pdom_e";
                }
                
                if (@{$lrange{$src}{$element}{$pdom_type}} > 1) {
                    {
                        no warnings; # Number::Range warns on EVERY single interger that overlaps
                        my $range = Number::Range->new(@{$lrange{$src}{$element}{$pdom_type}});
                        $union = $range->range;
                    }
                            
                    for my $r (split /\,/, $union) {
                        my ($ustart, $uend) = split /\.\./, $r;
                        my $seq = $self->subseq_pdoms($index, $src, $ustart, $uend);
                        my $k = join "_", $ustart, $uend;
                        $seqs{$k} = $seq;
                    }
                            
                    $self->concat_pdoms($src, $element, \%seqs, $fh);
                }
                else {
                    my ($nustart, $nuend, $str) = split /\|\|/, @{$pdoms{$src}{$element}{$pdom_type}}[0];
                    $self->subseq($index, $src, $element, $nustart, $nuend, $fh);
                }
                close $fh;
                %seqs   = ();
                %lrange = ();
                unlink $pdom_file if ! -s $pdom_file;
            }
        }
    }

    for my $file ($comp, $ppts, $pbs, $five_pr_ltrs, $three_pr_ltrs) {
        unlink $file if ! -s $file;
    }

    return $resdir
}

sub subseq_pdoms {
    my $self = shift;
    my ($index, $loc, $start, $end) = @_;

    my $location = "$loc:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);
    croak "\nERROR: Something went wrong. This is a bug, please report it.\n"
        unless $length;
    return $seq;
}

sub concat_pdoms {
    my $self = shift;
    my ($src, $elem, $seqs, $fh_out) = @_;
    my @ranges = map { split /\_/, $_ } keys %$seqs;
    my $start  = min(@ranges);
    my $end    = max(@ranges);
    my $id     = join "_", $elem, $src, $start, $end;

    my $concat_seq;
    for my $seq (values %$seqs) {
        $concat_seq .= $seq;
    }

    $concat_seq =~ s/.{60}\K/\n/g;
    say $fh_out join "\n", ">$id", $concat_seq;
}

sub collect_feature_args {
    my $self = shift;
    my ($dir) = @_;
    my (@fiveltrs, @threeltrs, @ppt, @pbs, @pdoms, %vmatch_args);
    find( sub { push @fiveltrs, $File::Find::name if -f and /5prime-ltrs.fasta$/ }, $dir);
    find( sub { push @threeltrs, $File::Find::name if -f and /3prime-ltrs.fasta$/ }, $dir);
    find( sub { push @ppt, $File::Find::name if -f and /ppts.fasta$/ }, $dir);
    find( sub { push @pbs, $File::Find::name if -f and /pbs.fasta$/ }, $dir);
    find( sub { push @pdoms, $File::Find::name if -f and /pdom.fasta$/ }, $dir);

    # ltr
    my $ltr5name = File::Spec->catfile($dir, 'dbcluster-5primeseqs');
    my $fiveargs = "-qspeedup 2 -dbcluster 80 0 $ltr5name -p -d -seedlength 30 ";
    $fiveargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 10000 -identity 80";
    $vmatch_args{fiveltr} = { seqs => \@fiveltrs, args => $fiveargs };

    my $ltr3name  = File::Spec->catfile($dir, 'dbcluster-3primeseqs');
    my $threeargs = "-qspeedup 2 -dbcluster 80 0 $ltr3name -p -d -seedlength 30 ";
    $threeargs .= "-exdrop 7 -l 80 -showdesc 0 -sort ld -best 10000 -identity 80";
    $vmatch_args{threeltr} = { seqs => \@threeltrs, args => $threeargs };

    # pbs/ppt
    my $pbsname = File::Spec->catfile($dir, 'dbcluster-pbs');
    my $pbsargs = "-dbcluster 90 90 $pbsname -p -d -seedlength 5 -exdrop 2 ";
    $pbsargs .= "-l 3 -showdesc 0 -sort ld -best 10000 -identity 90";
    $vmatch_args{pbs} = { seqs => \@pbs, args => $pbsargs, prefixlen => 1 };

    my $pptname = File::Spec->catfile($dir, 'dbcluster-ppt');
    my $pptargs = "-dbcluster 90 90 $pptname -p -d -seedlength 5 -exdrop 2 ";
    $pptargs .= "-l 3 -showdesc 0 -sort ld -best 10000 -identity 90";
    $vmatch_args{ppt} = { seqs => \@ppt, args => $pptargs, prefixlen => 5 };

    # pdoms
    my $pdomname = File::Spec->catfile($dir, 'dbcluster-pdoms');
    my $pdomargs = "-qspeedup 2 -dbcluster 80 0 $pdomname -p -d -seedlength 30 -exdrop 3 ";
    $pdomargs .= "-l 40 -showdesc 0 -sort ld -best 10000";
    $vmatch_args{pdoms} = { seqs => \@pdoms, args => $pdomargs };

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
    my $outfile = File::Spec->catfile($dir, 'all_vmatch_reports.txt');
    my $logfile = File::Spec->catfile($dir, 'all_vmatch_reports.log');
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $log, '>>', $logfile or die "\nERROR: Could not open file: $logfile\n";
    
    my $thr; 
    if ($threads % 3 == 0) {
	$thr = sprintf("%.0f", $threads/3);
    }
    elsif ($threads-1 % 3 == 0) {
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
    my $index = File::Spec->catfile($path, $name.'.index');
    my $vmrep = File::Spec->catfile($path, $name.'_vmatch-out.txt');
    my $log   = File::Spec->catfile($path, $name.'_vmatch-out.log');;

    my $mkvtreecmd = "mkvtree -db $db -dna -indexname $index -allout -v -pl ";
    if (defined $args->{$type}{prefixlen}) {
        $mkvtreecmd .= "$args->{$type}{prefixlen} ";
    }
    $mkvtreecmd .= "2>&1 > $log";
    my $vmatchcmd  = "vmatch $args->{$type}{args} $index > $vmrep";
    $self->run_cmd($mkvtreecmd);
    $self->run_cmd($vmatchcmd);
    unlink glob "$index*";
    unlink glob "$path/*.match";
    unlink $log;

    return $vmrep;
}

sub subseq {
    my $self = shift;
    my ($index, $loc, $elem, $start, $end, $out) = @_;

    my $location = "$loc:$start-$end";
    my ($seq, $length) = $index->get_sequence($location);
    croak "\nERROR: Something went wrong. This is a bug, please report it.\n"
        unless $length;

    my $id = join "_", $elem, $loc, $start, $end;

    $seq =~ s/.{60}\K/\n/g;
    say $out join "\n", ">$id", $seq;
}

sub parse_clusters {
    my $self = shift;
    my ($clsfile) = @_;
    my $genome = $self->genome;
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
        $name =~ s/$1//;
    }

    my (%cls, %all_seqs, %all_pdoms, $clusnum, $dom);
    open my $in, '<', $clsfile or die "\nERROR: Could not open file: $clsfile\n";

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
            $element =~ s/_\d+-?_?\d+$//;
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
                unlink $db;
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

    perldoc Tephra::Classify::LTRFams::Cluster


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
