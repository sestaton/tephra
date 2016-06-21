package Tephra::Genome::SoloLTRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use File::Copy;
use File::Path          qw(make_path remove_tree);
use IPC::System::Simple qw(capture system);
use Path::Class::File;
use Bio::DB::HTS::Kseq;
use Bio::AlignIO;
use Bio::SearchIO;
use Sort::Naturally;
use Set::IntervalTree;
use Carp 'croak';
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::Genome::SoloLTRSearch - Find solo-LTRs in a refence genome

=head1 VERSION

Version 0.03.1

=cut

our $VERSION = '0.03.1';
$VERSION = eval $VERSION;

has dir => (
      is       => 'ro',
      isa      => 'Path::Class::Dir',
      required => 1,
      coerce   => 1,
);

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has report => (
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

has seq => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

has percentident => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_percentident',
    lazy      => 1,
    default   => 0.39,
);

has percentcov => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_percentcov',
    lazy      => 1,
    default   => 0.80,
);

has matchlen => (
    is        => 'ro',
    isa       => 'Num',
    predicate => 'has_matchlen',
    lazy      => 1,
    default   => 80,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

sub find_soloLTRs {
    my $self = shift;
    my $dir = $self->dir;
    my $genome = $self->genome;
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $hmmsearch_summary = $self->report // File::Spec->catfile($path, $name.'_tephra_soloLTRs.tsv');
    my ($hmmbuild, $hmmsearch) = $self->_find_hmmer;

    print STDERR "Getting LTR alignments....";
    my $ltr_aln_files = $self->_get_ltr_alns($dir);
    say STDERR "done with alignments.";
    
    ## need masked genome here
    if (@$ltr_aln_files < 1 || ! -e $genome) {
	croak "\nERROR: No genome was found or the expected alignments files were not found. Exiting.";
    }
    
    print STDERR "Getting alignment statistics...";
    my $aln_stats = $self->_get_aln_len($ltr_aln_files); # return a hash-ref
    say STDERR "done with alignment statistics.";

    # make one directory
    my $model_dir = File::Spec->catdir($dir, 'Tephra_LTR_exemplar_models');
    unless ( -e $model_dir ) {
	make_path( $model_dir, {verbose => 0, mode => 0771,} );
    }
    
    print STDERR "Building LTR exemplar models...";
    for my $ltr_aln (@$ltr_aln_files) {
	$self->build_model($ltr_aln, $model_dir, $hmmbuild);	
	unlink $ltr_aln;	
    }
    say STDERR "done with models.";
    
    my @ltr_hmm_files;
    find( sub { push @ltr_hmm_files, $File::Find::name if -f and /\.hmm$/ }, $model_dir);
    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);

    print STDERR "Search genome with models...";
    for my $hmm (@ltr_hmm_files) {
	my $indiv_results_dir = $hmm;
	$indiv_results_dir =~ s/\.hmm.*$/\_search_out/;
	
	my $hmmsearch_out = $hmm;
	$hmmsearch_out =~ s/\.hmm.*$//;
	$hmmsearch_out .= '_'.$gname.'.hmmer';
	
	$self->search_with_models($hmmsearch_out, $hmm, $hmmsearch, $genome);
	    
	if (-e $hmmsearch_out) {   
	    $self->write_hmmsearch_report($aln_stats, $hmmsearch_out);
	} 
    }
    say STDERR "done searching with LTR models.";
    print STDERR "Collating results and writing GFF...";
    #my @zeroes;
    #find( sub { push @zeroes, $File::Find::name if -f and ! -s }, $model_dir );
    #unlink @zeroes;
    
    my @reports;
    find( sub { push @reports, $File::Find::name if -f and /\.txt$/ }, $model_dir );
    if (@reports) {
	$self->_collate(\@reports, $hmmsearch_summary);
    }
    else {
	say "\nWARNING: No solo-LTRs found so none will be reported. Exiting.\n";
	unlink $hmmsearch_summary if -e $hmmsearch_summary;
	exit(1);
    }

    $self->write_sololtr_gff($hmmsearch_summary);
    say STDERR "all done with solo-LTRs.";
}

sub write_sololtr_gff {
    my $self = shift;
    my ($hmmsearch_summary) = @_;
    my $genome  = $self->genome;
    my $outfile = $self->outfile;

    my $seqlen = $self->_get_seq_len($genome);
    open my $in, '<', $hmmsearch_summary or die "\nERROR: Could not open file: $hmmsearch_summary\n";
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";

    say $out "##gff-version 3";
    for my $id (nsort keys %$seqlen) {
	say $out join q{ }, "##sequence-region", $id, '1', $seqlen->{$id};
    }

    my $ct = 0;
    my (%features, %intervals, %ltr_lengths);
    while (my $line = <$in>) {
	chomp $line;
	next if $line =~ /^#/;
	my ($query, $query_length, $number_of_hits, $hit_name, $perc_coverage, $hsp_length, 
	    $hsp_perc_ident, $hsp_query_start, $hsp_query_end, $hsp_hit_start, $hsp_hit_end,
	    $search_type) = split /\t/, $line;
	$query =~ s/_ltrs_muscle-out//;

	if (exists $seqlen->{$hit_name}) {
	    $ct++;
	    say $out join "\t", $hit_name, 'Tephra', 'solo_LTR', $hsp_hit_start, $hsp_hit_end, 
	        '.', '?', '.', "ID=solo_LTR$ct;Parent=$query;Name=solo_LTR;Ontology_term=SO:0001003";
	}
	else {
	    die "\nERROR: $hit_name not found in $genome. This should not happen. Exiting.\n";
	}
    }
    close $in;
    close $out; 
}

sub build_model {
    my $self = shift;
    my ($aln, $model_dir, $hmmbuild) = @_;

    my $hmmname = basename($aln);
    $hmmname =~ s/\.aln$/\.hmm/;
    my $model_path = File::Spec->catfile($model_dir, $hmmname);
    # hmmer3 has --dna and --cpu options and not -g (global)
    my $hmm_cmd = "$hmmbuild -f --nucleic $model_path $aln";

    $self->capture_cmd($hmm_cmd);
}

sub search_with_models {
    my $self = shift;
    my ($search_out, $hmm_model, $hmmsearch, $fasta) = @_;

    ## no -o option in hmmer2
    my $hmmsearch_cmd = "$hmmsearch --cpu 12 $hmm_model $fasta > $search_out";
    $self->run_cmd($hmmsearch_cmd);
    unlink $hmm_model; # if $self->clean;
}

sub _get_ltr_alns {
    my $self = shift;
    my ($dir) = @_;
    my (@ltrseqs, @aligns);

    find( sub { push @ltrseqs, $File::Find::name if -f and /exemplar_ltrs.fasta$/ }, $dir);

    for my $ltrseq (@ltrseqs) {
	my ($name, $path, $suffix) = fileparse($ltrseq, qr/\.[^.]*/);
	my $tre = File::Spec->catfile($path, $name.'.dnd');
	my $aln = File::Spec->catfile($path, $name.'_muscle-out.aln');
	my $log = File::Spec->catfile($path, $name.'_muscle-out.log');
	
	my $clwcmd = "muscle -clwstrict -in $ltrseq -out $aln 2>$log";
	$self->capture_cmd($clwcmd);
	unlink $log;
	push @aligns, $aln;
    }

    return \@aligns;
}

sub _collate {
    my $self = shift;
    my ($files, $outfile) = @_;

    ##TODO: use lexical vars instead for array indexes
    my (%seen, %parsed_alns);
    for my $file (@$files) {
	open my $fh_in, '<', $file or die "\nERROR: Could not open file: $file\n";
	while (my $line = <$fh_in>) {
	    chomp $line;
	    next if $line =~ /^#/;
	    my @f = split /\t/, $line;
	    if (exists $seen{$f[3]}{$f[9]}) {
		my ($len, $pid) = split /\|\|/, $seen{$f[3]}{$f[9]};
		if ($f[10]-$f[9] > $len && $f[6] > $pid) {
		    $parsed_alns{$f[3]}{$f[9]} =  join "||", @f;
		}
	    }
	    $parsed_alns{$f[3]}{$f[9]} =  join "||", @f;
	    $seen{$f[3]}{$f[9]} = join "||", $f[10]-$f[9], $f[6];
	}
	close $fh_in;
    }

    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    say $out join "\t", "#query", "query_length", "number_of_hits", "hit_name",
        "perc_coverage","hsp_length", "hsp_perc_ident","hsp_query_start", "hsp_query_end",
        "hsp_hit_start", "hsp_hit_end","search_type";

    for my $chr (nsort keys %parsed_alns) {
	for my $start (sort { $a <=> $b } keys %{$parsed_alns{$chr}}) {
	    my @f = split /\|\|/, $parsed_alns{$chr}{$start};
	    say $out join "\t", @f;
	}
    }
    close $out;
}

sub _find_hmmer {
    my $self = shift; 

    my $config      = Tephra::Config::Exe->new->get_config_paths;
    my ($hmmer2bin) = @{$config}{qw(hmmer2bin)};
    my $hmmbuild    = File::Spec->catfile($hmmer2bin, 'hmmbuild');
    my $hmmsearch   = File::Spec->catfile($hmmer2bin, 'hmmsearch');

    if (-e $hmmbuild && -x $hmmbuild &&
	-e $hmmsearch && -x $hmmsearch) {
	return ($hmmbuild, $hmmsearch);
    }
    else {
	croak "\nERROR: Could not get HMMERv2 PATH. This indicates that Tephra was not configured correctly. Exiting.\n";
    }
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

sub _get_aln_len {
    my $self = shift;
    my ($aln_files) = @_;
    
    my %aln_stats;

    for my $aligned (@$aln_files) {
	my $ltr_retro = basename($aligned);
	$ltr_retro =~ s/\.aln//;
	my $aln_in = Bio::AlignIO->new(-file => $aligned, -format => 'clustalw');

	while ( my $aln = $aln_in->next_aln() ) {
	    $aln_stats{$ltr_retro} = $aln->length;
	}	
    }

    return \%aln_stats;
}

sub write_hmmsearch_report {
    my $self = shift;
    my $genome     = $self->genome;
    my $match_pid  = $self->percentident;
    my $match_len  = $self->matchlen;
    my $match_pcov = $self->percentcov;
    my ($aln_stats, $search_report) = @_;

    my $parsed = $search_report;
    $parsed =~ s/\.hmmer$/\_hmmer_parsed.txt/;
    my $seqfile = $parsed;
    $seqfile =~ s/\.txt$/\_seq.fasta/;
    my $element = $search_report;

    my $model_type = 'local';

    my ($gname, $gpath, $gsuffix) = fileparse($genome, qr/\.[^.]*/);
    my ($name, $path, $suffix)    = fileparse($search_report, qr/\.[^.]*/);
    $element =~ s/_$gname*//;
    open my $out, ">$parsed" or die "\nERROR: Could not open file: $parsed\n";

    my $seq;
    if (defined $self->seq) {
	open $seq, ">$seqfile" or die "\nERROR: Could not open file: $seqfile";
    }

    my $hmmerin = Bio::SearchIO->new(-file => $search_report, -format => 'hmmer');

    my ($positions, $matches) = (0, 0);
    while ( my $result = $hmmerin->next_result() ) {
	my $query    = $result->query_name();
	my $num_hits = $result->num_hits();
	my $qlen     = $aln_stats->{$query};
	die "\nERROR: Could not determine query length for: $query" 
	    unless defined $qlen;
	while ( my $hit = $result->next_hit() ) {
	    my $hitid = $hit->name();
	    while ( my $hsp = $hit->next_hsp() ) {
		my $percent_q_coverage = sprintf("%.2f", $hsp->length('query')/$qlen * 100);
		my @ident_pos = $hsp->seq_inds('query','identical');
		my $hspgaps   = $hsp->gaps;
		my $hsplen    = $hsp->length('total');
		my $hstart    = $hsp->start('hit');
		my $hstop     = $hsp->end('hit');
		my $qstart    = $hsp->start('query');
		my $qstop     = $hsp->end('query');
		my $qstring   = $hsp->query_string;
		my $qhsplen   = $hsp->length('query');
		
		if (exists $aln_stats->{$query}) {
		    my $percent_coverage = sprintf("%.2f",$hsplen/$aln_stats->{$query});
		    $positions++ for @ident_pos;
		    my $pid = sprintf("%.2f", $positions/$aln_stats->{$query} * 100);
		    if ( $hsplen >= $match_len && $hsplen >= $aln_stats->{$query} * $match_pcov ) {
			if ($pid >= $match_pid) {
			    $matches++;
			    say $out join "\t", $query, $aln_stats->{$query}, $num_hits, $hitid, 
			        $percent_q_coverage, $hsplen, $pid, $qstart, $qstop, $hstart, 
			        $hstop, $model_type;

			    if ($self->seq) {
				my $seqid = '>'.$query.'_'.$hitid.'_'.$hstart.'_'.$hstop; 
				## It makes more sense to show the location of the hit
				## Also, this would pave the way for creating a gff of solo-LTRs
				## my $seqid = ">".$query."|".$hitid."_".$hstart."-".$hstop
				say $seq join "\n", $seqid, $qstring;
			    }
			}		
		    }
		    $positions = 0;
		}
	    }
	}
    }
    close $out;
    close $seq;
    unlink $seqfile, $parsed if $matches == 0;

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
