package Tephra::Genome::FragmentSearch;

use 5.014;
use Moose;
use File::Spec;
use File::Basename;
use Cwd                 qw(abs_path);
use IPC::System::Simple qw(capture system);
use Time::HiRes         qw(gettimeofday);
use Sort::Naturally;
use Set::IntervalTree;
use Parallel::ForkManager;
use Carp 'croak';
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::Util',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Genome::FragmentSearch - Find fragmented transposons in a refence genome

=head1 VERSION

Version 0.09.2

=cut

our $VERSION = '0.09.2';
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
    default   => 200,
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

    my %reports;
    my $t0 = gettimeofday();
    my $blast_report = $self->search_genome($threads);
    my $blast_files  = $self->split_refs($blast_report);
    my $seqlen       = $self->_get_seq_len($genome);

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my $logfile = File::Spec->catfile($gpath, 'tephra_fragment_searches.log');
    open my $log, '>>', $logfile or die "\nERROR: Could not open file: $logfile\n";
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";

    if (! -s $outfile) {
	say $out '##gff-version 3';
	for my $id (nsort keys %$seqlen) {
	    say $out join q{ }, '##sequence-region', $id, '1', $seqlen->{$id};
	}
    }

    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        warn "Caught SIGINT; Waiting for child processes to finish.";
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      for my $src (nsort keys %$data_ref) {
				  my $windows = $data_ref->{$src};
				  $self->write_fragment_gff($src, $windows, $out);
				  undef $windows;
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

    unlink $blast_report;
    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $log "\n========> Finished searching ",scalar(@$blast_files)," sequence in $genome in $final_time minutes.";
    close $log;
}

sub collapse_overlaps {
    my $self = shift;
    my ($file) = @_;

    my $tree = Set::IntervalTree->new;
    my %windows;

    open my $l, '<', $file or die "\nERROR: Could not open file: $file\n";

    my $line = <$l>;
    chomp $line;
    my @f = split /\t/, $line;

    my %trees;
    $tree->insert({ match => $f[0], start => $f[8], end => $f[9], len => $f[3] }, 
                  $f[8], $f[9]);
    
    $windows{$f[8]} =
        { match => $f[0], start => $f[8], end => $f[9], len => $f[3], evalue => $f[10], strand => $f[12] };
    
    while (my $line = <$l>) {
        chomp $line;
        
        my ($queryId, $subjectId, $percIdentity, $alnLength, $mismatchCount, $gapOpenCount, $queryStart, 
            $queryEnd, $subjectStart, $subjectEnd, $eVal, $bitScore, $strand) = split /\t/, $line;
        
        my $res = $tree->fetch($subjectStart, $subjectEnd);

        if (@$res) {
            #dd $res and exit
                #if @$res > 1;
            for my $overlap (@$res) {
                my ($ostart, $oend, $match, $subj, $olen) = @{$overlap}{qw(start end match id len)};
                my $overl = $subjectEnd - $oend;
                next if $subjectStart >= $ostart && $subjectEnd <= $oend;
                my $nlen = $subjectEnd - $ostart + 1;
                
                $windows{$ostart} =
                    { match => $queryId, start => $ostart, end => $subjectEnd, len => $nlen, evalue => $eVal, strand => $strand};
                
                $tree->remove($ostart, $oend);
                $tree->insert({ id => $subjectId, match => $queryId, start => $ostart, end => $subjectEnd, len => $nlen }, 
                              $ostart, $subjectEnd);
                #$tree->remove($ostart, $oend); 
            }
        }
	else {
            $tree->insert({ match => $queryId, start => $subjectStart, end => $subjectEnd, len => $alnLength }, 
                          $subjectStart, $subjectEnd);
            
            $windows{$subjectStart} =
                { match => $queryId, start => $subjectStart, end => $subjectEnd, len => $alnLength, 
                  evalue => $eVal, strand => $strand };
        }
    }
    close $l;

    undef $tree;
    my $src = $file;
    $src =~ s/_tmp.*//;
    unlink $file;

    return ($src, \%windows);
}

sub write_fragment_gff {
    my $self = shift;
    my ($src, $windows, $out) = @_;

    my $filt = 0;
    for my $s (sort { $a <=> $b } keys %$windows) {
        $filt++;
        my ($match, $sstart, $send, $slen, $seval, $sstr) = @{$windows->{$s}}{qw(match start end len evalue strand)};
        say $out join "\t", $src, 'BLASTN', 'similarity', $sstart, $send, $seval, $sstr, '.', "ID=$match"."_fragment_$src-$filt"; 
    }

    return;
}

sub search_genome {
    my $self = shift;
    my ($thr) = @_;
    my $genome   = $self->genome->absolute->resolve;
    my $repeatdb = $self->repeatdb->absolute->resolve;

    my ($dbname, $dbpath, $dbsuffix) = fileparse($repeatdb, qr/\.[^.]*/);
    my ($qname, $qpath, $qsuffix)    = fileparse($genome, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile( abs_path($qpath), $qname."_$dbname".'.bln' );

    my $blastdb = $self->make_blastdb($genome);
    my $blast_report = $self->run_blast({ query => $repeatdb, db => $blastdb, outfile => $outfile, 
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

    open my $b, '<', $blast or die "\nERROR: Could not open file: $blast\n";

    my (%refs, %coords, @outfiles);

    my (%ofhs); #, $ofh);
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
            #open my $ofh, '|-', "sort -nk9,9 >$outf" or die $!;
	    open my $ofh, '>', $outf or die "\nERROR: Could not open file: $outf\n";
            push @outfiles, $outf;
            $ofhs{$f[1]} = $ofh;
            say $ofh join "\t", @f;
        }
    }
    close $_ for keys %ofhs;

    return \@outfiles;
}

sub _get_seq_len {
    my $self = shift;
    my ($genome) = @_;
    
    my %len;

    my $kseq = Bio::DB::HTS::Kseq->new($genome);
    my $iter = $kseq->iterator();

    while ( my $seq = $iter->next_seq() ) {
	my $id  = $seq->name;
	my $seq = $seq->seq;
	$len{$id} = length($seq);
    }       

    return \%len;
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

S. Evan Staton, C<< <statonse at gmail.com> >>

=head1 BUGS

Please report any bugs or feature requests through the project site at 
L<https://github.com/sestaton/tephra/issues>. I will be notified,
and there will be a record of the issue. Alternatively, I can also be 
reached at the email address listed above to resolve any questions.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Tephra::Genome::SoloLTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
