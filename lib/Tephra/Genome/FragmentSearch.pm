package Tephra::Genome::FragmentSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Basename;
use File::Temp          qw(tempfile);
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(capture system);
use Time::HiRes         qw(gettimeofday);
use Sort::Naturally;
use Parallel::ForkManager;
use Carp 'croak';
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::File',
     'Tephra::Role::Util',
     'Tephra::Role::Run::Any',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Genome::FragmentSearch - Find fragmented transposons in a refence genome

=head1 VERSION

Version 0.14.0

=cut

our $VERSION = '0.14.0';
$VERSION = eval $VERSION;

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has repeatdb => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 0,
    coerce   => 1,
);

has outfile => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has percentid => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_percentid',
    lazy      => 1,
    default   => 80,
);

has matchlen => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_matchlen',
    lazy      => 1,
    default   => 100,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

sub find_transposon_fragments {
    my $self = shift;
    my $threads  = $self->threads;
    my $genome   = $self->genome->absolute->resolve;
    my $repeatdb = $self->repeatdb->absolute->resolve;
    my $outfile  = $self->outfile; 

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my $logfile = File::Spec->catfile($gpath, 'tephra_fragment_searches.log');

    open my $log, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";

    my (%reports, %window_refs);
    my $t0 = gettimeofday();
    my $index = $self->index_ref($genome);

    my ($seqlen, $genome_parts) = $self->split_genome($genome);
    my ($rdb, $rdb_is_compressed) = $self->_make_temp_repeatdb;

    for my $part (@$genome_parts) { 
	#say STDERR join q{ }, "\npart -> $part ; repeatdb -> $rdb";
	my $blastdb = $self->make_blastdb($part);
	#say STDERR "\nBLASTDB $blastdb";
	my $blast_report = $self->search_genome($rdb, $blastdb, $threads);
	#say STDERR "\nBLASTREP $blast_report";

	my $blast_files = $self->split_refs($blast_report);
	
	my $pm = Parallel::ForkManager->new($threads);
	local $SIG{INT} = sub {
	    warn "Caught SIGINT; Waiting for child processes to finish.";
	    $pm->wait_all_children;
	    exit 1;
	};
	
	$pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
				  for my $src (nsort keys %$data_ref) {
				      $window_refs{$src} = $data_ref->{$src};
				  }
				  my $t1 = gettimeofday();
				  my $elapsed = $t1 - $t0;
				  my $time = sprintf("%.2f",$elapsed/60);
				  say $log basename($ident),
				  " just finished with PID $pid and exit code: $exit_code in $time minutes";
			    } );
	
	for my $file (nsort @$blast_files) {
	    $pm->start($file) and next;
	    $SIG{INT} = sub { $pm->finish };
	    
	    my ($src, $windows) = $self->collapse_overlaps($file);
	    $reports{$src} = $windows;
	    
	    $pm->finish(0, \%reports);
	}
	
	$pm->wait_all_children;
	unlink $blastdb, $blast_report;
    }

    ## write output
    my ($oname, $opath, $osuffix) = fileparse($outfile, qr/\.[^.]*/);
    my $fafile  = File::Spec->catfile($opath, $oname.'.fasta');

    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    open my $faout, '>>', $fafile or die "\n[ERROR]: Could not open file: $fafile\n";

    if (! -s $outfile) {
        say $out '##gff-version 3';
        for my $id (nsort keys %$seqlen) {
            say $out join q{ }, '##sequence-region', $id, '1', $seqlen->{$id};
        }
    }

    for my $src (nsort keys %window_refs) {
	$self->write_fragment_gff($index, $src, $window_refs{$src}, $out, $faout);
    }

    close $out;
    close $faout;

    ## clean up
    my $exe_conf = Tephra::Config::Exe->new->get_config_paths;
    my ($gt, $cleanpp) = @{$exe_conf}{qw(gt cleanpp)};

    $self->capture_cmd($cleanpp);
    $self->capture_cmd($gt, 'clean');
    unlink $genome.'.fai';
    unlink $_ for @$genome_parts;
    unlink $rdb if $rdb_is_compressed;

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished searching ",scalar(keys %$seqlen)," sequence in $genome in $final_time minutes.";
    close $log;

    return;
}

sub collapse_overlaps {
    my $self = shift;
    my ($file) = @_;
    my $matchlen = $self->matchlen;

    my %windows;

    open my $l, '<', $file or die "\n[ERROR]: Could not open file: $file\n";
    
    my (@f, $prev_start, $prev_end, $prev_strand, $prev_len);
    line : { 
	my $line = <$l>;
	chomp $line;
	@f = split /\t/, $line;

	redo line unless $f[3] >= $matchlen;
	
	$prev_len = $f[9] - $f[8] + 1;
	$prev_strand = $f[12];
	($prev_start, $prev_end) = @f[8,9];
    }
    
    redo line unless defined $prev_start && defined $prev_end;
    $windows{$f[8]} =
        { match => $f[0], start => $prev_start, end => $prev_end, len => $prev_len, evalue => $f[10], strand => $prev_strand };

    while (my $line = <$l>) {
        chomp $line;
        
        my ($queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, 
            $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore, $strand) = split /\t/, $line;
     
	my $aln_length = $subjectEnd - $subjectStart + 1;
	next unless $aln_length >= $matchlen;

	if ($subjectStart > $prev_end) { 
	    $windows{$subjectStart} =
	        { match => $queryId, start => $subjectStart, end => $subjectEnd, len => $aln_length, evalue => $eVal, strand => $strand };

            ($prev_start, $prev_end) = ($subjectStart, $subjectEnd);
        }
    }
    close $l;

    my $src = $file;
    $src =~ s/_tmp.*//;
    unlink $file;
    
    return ($src, \%windows);
}

sub write_fragment_gff {
    my $self = shift;
    my ($index, $src, $windows, $out, $faout) = @_;
    my $matchlen = $self->matchlen;

    my $filt = 0;
    for my $s (sort { $a <=> $b } keys %$windows) {
        $filt++;
        my ($match, $sstart, $send, $slen, $seval, $sstr) = @{$windows->{$s}}{qw(match start end len evalue strand)};
	my ($seq, $length) = $self->get_full_seq($index, $src, $sstart, $send);
	$seq =~ s/.{60}\K/\n/g;
	my $id = join "_", $match, 'fragment', "$src-$filt";
	my $seqid = join "_", $id, $src, $sstart, $send;

	say $faout join "\n", ">".$seqid, $seq;
	say $out join "\t", $src, 'BLASTN', 'similarity', $sstart, $send, $seval, $sstr, '.', "ID=$id";
    }

    return;
}

sub split_genome {
    my $self = shift; 
    my ($genome) = @_;

    my $numreads = 500;
    my $num    = 0;
    my $count  = 0;
    my $fcount = 1;

    my ($iname, $ipath, $isuffix) = fileparse($genome, qr/\.[^.]*/);
    if ($iname =~ /\.fa.*/) {
	($iname, $ipath, $isuffix) = fileparse($iname, qr/\.[^.]*/);
    }

    my ($out, @split_files, %len);
    
    my $fname = File::Spec->catfile($ipath, $iname."_".$fcount.$isuffix);
    open $out, '>', $fname or die "\n[ERROR]: Could not open file: $fname\n";
    
    push @split_files, $fname;

    my $kseq = Bio::DB::HTS::Kseq->new($genome);
    my $iter = $kseq->iterator();

    while ( my $seqobj = $iter->next_seq() ) {
        my $id  = $seqobj->name;
        my $seq = $seqobj->seq;
	$len{$id} = length($seq);

	if ($count % $numreads == 0 && $count > 0) {
	    $fcount++;
	    #$fname = $iname."_".$fcount.$isuffix;
	    $fname = File::Spec->catfile($ipath, $iname."_".$fcount.$isuffix);
	    open $out, '>', $fname or die "\n[ERROR]: Could not open file: $fname\n";
	    
	    push @split_files, $fname;
	}
	else {
	    say $out join "\n", ">".$id, $seq;
	}
	$count++;
    }
    close $out;

    return (\%len, \@split_files);
}

sub search_genome {
    my $self = shift;
    my ($rdb, $blastdb, $thr) = @_;
    my $genome   = $self->genome->absolute->resolve;
    #my $repeatdb = $self->repeatdb->absolute->resolve;

    my ($dbname, $dbpath, $dbsuffix) = fileparse($rdb, qr/\.[^.]*/);
    my ($qname, $qpath, $qsuffix)    = fileparse($genome, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile( abs_path($qpath), $qname."_$dbname".'.bln' );

    my $blast_report = $self->run_blast({ query => $rdb, db => $blastdb, outfile => $outfile, 
					  threads => $thr, evalue => 1e-10, sort => 'coordinate' });
    my @dbfiles = glob "$blastdb*";
    unlink @dbfiles;

    return $blast_report;
}

sub split_refs {
    my $self = shift;
    my ($blast) = @_;
    my $matchlen  = $self->matchlen;
    my $percentid = $self->percentid;

    open my $b, '<', $blast or die "\n[ERROR]: Could not open file: $blast\n";

    my (%refs, %coords, @outfiles, %ofhs);
    while (my $line = <$b>) {
        chomp $line;
        my @f = $self->_format_hits($line);
        next unless @f; # skip hit if not properly formatted
        next unless $f[2] >= $percentid && $f[3] >= $matchlen; 

        if (exists $ofhs{$f[1]}) {
            my $ofh = $ofhs{$f[1]};
            say $ofh join "\t", @f;
        }
        else {
            my $outf = $f[1].'_tmp.bln';
            # hits have to be sorted by reference for interval tree to work
            open my $ofh, '|-', "sort -nk9,9 >$outf" or die $!;
	    #open my $ofh, '>', $outf or die "\n[ERROR]: Could not open file: $outf\n";
            push @outfiles, $outf;
            $ofhs{$f[1]} = $ofh;
            say $ofh join "\t", @f;
        }
    }
    close $_ for keys %ofhs;

    return \@outfiles;
}

sub _make_temp_repeatdb {
    my $self = shift;
    my $rdb = $self->repeatdb->absolute->resolve;

    my $rdb_is_compressed = 0;
    if ($rdb =~/\.gz$|\.bz2$/) {
        my $fh = $self->get_fh($rdb);

        #$fasta =~ s/\.gz$|\.bz2$//;
        $rdb =~ s/\.gz$|\.bz2$//;
        my ($dbname, $dbpath, $dbsuffix) = fileparse($rdb, qr/\.[^.]*/);
        #my ($faname, $fapath, $fasuffix) = fileparse($fasta, qr/\.[^.]*/);
        #$outfile = File::Spec->catfile($fapath, $faname.'_'.$dbname.'.bln');
        
        my $tmpiname = $dbname.'_XXXX';
        my ($tmp_fh, $tmp_rdb) = tempfile( TEMPLATE => $tmpiname, DIR => $dbpath, SUFFIX => '', UNLINK => 0 );

        while (my $line = <$fh>) {
            chomp $line;
            say $tmp_fh $line;
        }
        close $fh;
        close $tmp_fh;

        $rdb = $tmp_rdb;
        #$blastdb = $self->make_blastdb($tmp_rdb);
	$rdb_is_compressed = 1;
    }
    #else {
        #my ($dbname, $dbpath, $dbsuffix) = fileparse($rdb, qr/\.[^.]*/);
        #my ($faname, $fapath, $fasuffix) = fileparse($fasta, qr/\.[^.]*/);

        #$outfile = File::Spec->catfile($fapath, $faname.'_'.$dbname.'.bln');
    #}

    return ($rdb, $rdb_is_compressed);
}

sub _format_hits {
    my $self = shift;
    my ($line) = @_;

    my ($start, $stop, $strand);
    my @f = split /\t/, $line;
    return 0 unless @f == 12;

    if ($f[8] > $f[9]) {
        $start = $f[9];
        $stop = $f[8];
        $strand = '-';
    }
    else {
        $start = $f[8];
        $stop = $f[9];
        $strand = '+';
    }

    @f = (@f[0..7], $start, $stop, @f[10..11], $strand);

    return @f;
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

    perldoc Tephra::Genome::FragmentSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
