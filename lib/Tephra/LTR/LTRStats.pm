package Tephra::LTR::LTRStats;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use Sort::Naturally;
use File::Spec;
use File::Find;
use File::Basename;
use File::Path      qw(make_path remove_tree);
use File::Copy      qw(move copy);
use List::MoreUtils qw(indexes any);
use Time::HiRes     qw(gettimeofday);
use Log::Any        qw($log);
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use Bio::Tools::GFF;
use Parallel::ForkManager;
use Cwd;
use Try::Tiny;
use Tephra::Config::Exe;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::PAML';

=head1 NAME

Tephra::LTR::LTRStats - Calculate the age distribution of LTR retrotransposons

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

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

#
# methods
#
sub extract_ltr_features {
    my $self = shift;
    my $fasta = $self->genome;
    my $gff   = $self->gff;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $dir = File::Spec->catdir($name.'_ltrages');
    unless ( -d $dir ) {
	make_path( $dir, {verbose => 0, mode => 0771,} );
    }

    my $gffio = Bio::Tools::GFF->new( -file => $gff, -gff_version => 3 );

    my ($start, $end, $elem_id, $family, $key, %feature, %ltrs, %seen);
    while (my $feature = $gffio->next_feature()) {
	if ($feature->primary_tag eq 'LTR_retrotransposon') {
	    my @string = split /\t/, $feature->gff_string;
	    ($family, $elem_id) = ($string[8] =~ /Family=?\s+?(RL._.*?family\d+)\s+\;\s+ID=?\s+?(LTR_retrotransposon\d+)\s+?/);
	    ($start, $end) = ($feature->start, $feature->end);
	    $key = join "||", $family, $elem_id, $start, $end;
	}
	next unless defined $start && defined $end;
	if ($feature->primary_tag eq 'long_terminal_repeat') {
	    my @string = split /\t/, $feature->gff_string;
	    if ($feature->start >= $start && $feature->end <= $end) {
		my $ltrkey = join "||", $string[0], $feature->primary_tag, @string[3,4,6];
		push @{$ltrs{$key}{'ltrs'}}, $ltrkey unless exists $seen{$ltrkey};
		$seen{$ltrkey} = 1;
	    }
	}
    }

    my %pdoms;
    my $ltrct = 0;
    for my $ltr (sort keys %ltrs) {
	my ($family, $element, $rstart, $rend) = split /\|\|/, $ltr;
	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
	    my $ltr_file = join "_", $family, $element, $src, $rstart, $rend, 'ltrs.fasta';
	    my $ltrs_out = File::Spec->catfile($dir, $ltr_file);
	    open my $ltrs_outfh, '>>', $ltrs_out or die "\nERROR: Could not open file: $ltrs_out\n";
	    my $lfname = $element;
	    my $orientation;
	    if ($ltrct) {
		$orientation = '5prime' if $strand eq '+';
		$orientation = '3prime' if $strand eq '-';
		$orientation = 'unknown-l' if $strand eq '.';
		$lfname .= "_$orientation-ltr.fasta";
		my $fiveprime_tmp = File::Spec->catfile($dir, $lfname);
		$self->subseq($fasta, $src, $family, $element, $s, $e, $fiveprime_tmp, $ltrs_outfh, $orientation);
		$ltrct = 0;
	    }
	    else {
		$orientation = '3prime' if $strand eq '+';
		$orientation = '5prime' if $strand eq '-';
		$orientation = 'unknown-r' if $strand eq '.';
		$lfname .= "_$orientation-ltr.fasta";
		my $threeprime_tmp = File::Spec->catfile($dir, $lfname);
		$self->subseq($fasta, $src, $family, $element, $s, $e, $threeprime_tmp, $ltrs_outfh, $orientation);
		$ltrct++;
	    }
	    close $ltrs_outfh;
	}
    }

    return $dir;
}

sub collect_feature_args {
    my $self = shift;
    my ($dir) = @_;

    my (@ltrs, %aln_args);
    my $wanted  = sub { push @ltrs, $File::Find::name if -f && /ltrs.fasta$/ };
    my $process = sub { grep ! -d, @_ };
    find({ wanted => $wanted, preprocess => $process }, $dir);
    
    $aln_args{ltrs} = { seqs => \@ltrs };

    return \%aln_args;
}

sub align_features {
    my $self = shift;
    my ($dir) = @_;
    my $threads = $self->threads;
    my $outfile = $self->outfile;

    my $resdir = File::Spec->catdir($dir, 'divergence_time_stats');
    
    unless ( -d $resdir ) {
	make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $args = $self->collect_feature_args($dir);
    #dd $args; ## debug

    my $t0 = gettimeofday();
    my $ltrrts = 0;
    my $logfile = File::Spec->catfile($dir, 'all_aln_reports.log');
    open my $logfh, '>>', $logfile or die "\nERROR: Could not open file: $logfile\n";
    
    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        $log->warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      my $file = $$data_ref;
			      my $t1 = gettimeofday();
			      my $elapsed = $t1 - $t0;
			      my $time = sprintf("%.2f",$elapsed/60);
			      say $logfh basename($ident)," just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    for my $type (keys %$args) {
	for my $db (@{$args->{$type}{seqs}}) {
	    $ltrrts++;
	    $pm->start($db) and next;
	    $SIG{INT} = sub { $pm->finish };
	    $self->process_align_args($db, $resdir);

	    $pm->finish(0, \$db);
	}
    }
    $pm->wait_all_children;

    my @agefiles;
    find( sub { push @agefiles, $File::Find::name if -f and /divergence.txt$/ and -s }, $resdir );

    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    say $out join "\t", "LTR-ID", "Divergence", "Age", "Ts:Tv";

    for my $file (@agefiles) {
	$self->collate($file, $out);
    }
    close $out;

    if ($self->clean) {
	my @alnfiles;
	my $wanted  = sub { push @alnfiles, $File::Find::name 
				if -f && /ltrs.fasta$|\.aln$|\.log$|\.phy$|\.dnd$|\.ctl$|\-divergence.txt$/ };
	my $process = sub { grep ! -d, @_ };
	find({ wanted => $wanted, preprocess => $process }, $dir);
	unlink @alnfiles;

	remove_tree( $dir, { safe => 1 } );
    }

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $logfh "\n========> Finished calculating insertion time for $ltrrts LTR-RTs in $final_time minutes";
    close $logfh;

    return;
}

sub process_baseml_args {
    my $self = shift;
    my ($phy, $dnd, $resdir) = @_;

    my $cwd = getcwd();
    my ($pname, $ppath, $psuffix) = fileparse($phy, qr/\.[^.]*/);
    my $divfile = File::Spec->catfile($ppath, $pname.'-divergence.txt');
    $divfile = basename($divfile);

    my $divergence = $self->_check_divergence($phy);

    if ($divergence > 0) {
	my $baseml_args = $self->create_baseml_files({ phylip => $phy, treefile => $dnd });
	$self->run_baseml({ working_dir => $ppath, current_dir => $cwd }); 
	$self->parse_baseml({ divergence_file => $divfile,
			      results_dir     => $resdir,
			      treefile        => $dnd,
			      phylip          => $phy,
			      outfile         => $baseml_args->{outfile},
			      control_file    => $baseml_args->{control_file} });
    }
    else {
	my $element = basename($phy);
	$element =~ s/_ltrs_clustal-out.*//;
	open my $divout, '>', $divfile or die "\nERROR: Could not open divergence file: $divfile\n";
	say $divout join "\t", $element, $divergence , '0', '0';
	close $divout;
	my $dest_file = File::Spec->catfile($resdir, $divfile);
	copy($divfile, $dest_file) or die "\nERROR: Copy failed: $!";
	unlink $divfile;
    }
}

sub process_align_args {
    my $self = shift;
    my ($db, $resdir) = @_;

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $pdir = File::Spec->catdir($path, $name.'_pamltmp');
    make_path( $pdir, {verbose => 0, mode => 0771,} );

    my $fas = File::Spec->catfile($pdir, $name.$suffix);
    copy($db, $fas) or die "\nERROR: Copy failed: $!";
    my $tre = File::Spec->catfile($pdir, $name.'.dnd');
    my $aln = File::Spec->catfile($pdir, $name.'_clustal-out.aln');
    my $dnd = File::Spec->catfile($pdir, $name.'_clustal-out.dnd');
    my $log = File::Spec->catfile($pdir, $name.'_clustal-out.log');

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($clustalw2) = @{$config}{qw(clustalw)};
    my $clwcmd = "$clustalw2 -infile=$fas -outfile=$aln 2>$log";
    $self->capture_cmd($clwcmd);
    my $phy = $self->parse_aln($aln, $tre, $dnd);
    $self->process_baseml_args($phy, $dnd, $resdir);
    
    return;
}

sub parse_aln {
    my $self = shift;
    my ($aln, $tre, $dnd) = @_;

    my ($name, $path, $suffix) = fileparse($aln, qr/\.[^.]*/);
    my $phy = File::Spec->catfile($path, $name.'.phy');
    
    my $aln_in  = Bio::AlignIO->new(-file  => $aln,    -format => 'clustalw');
    my $aln_out = Bio::AlignIO->new(-file  => ">$phy", -format => 'phylip', -flag_SI => 1, -idlength => 20);

    while (my $alnobj = $aln_in->next_aln) {
	$aln_out->write_aln($alnobj);
    }

    my $tre_in  = Bio::TreeIO->new(-file => $tre,    -format => 'newick');
    my $tre_out = Bio::TreeIO->new(-file => ">$dnd", -format => 'newick');

    while (my $treobj = $tre_in->next_tree) {
	for my $node ($treobj->get_nodes) {
	    my $id = $node->id;
	    next unless defined $id;
	    my $newid = substr $id, 0, 20;
	    $node->id($newid);
	}
	$tre_out->write_tree($treobj);
    }
    unlink $tre;
    
    return $phy;
}

sub subseq {
    my $self = shift;
    my ($fasta, $loc, $family, $elem, $start, $end, $tmp, $out, $orient) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($samtools) = @{$config}{qw(samtools)};
    my $cmd = "$samtools faidx $fasta $loc:$start-$end > $tmp";
    $self->run_cmd($cmd);

    my $id = join "_", $orient, $loc, $elem, "$start-$end";
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

sub collate {
    my $self = shift;
    my ($file_in, $fh_out) = @_;
    my $lines = do { 
	local $/ = undef; 
	open my $fh_in, '<', $file_in or die "\nERROR: Could not open file: $file_in\n";
	<$fh_in>;
    };
    print $fh_out $lines;
}

sub _check_divergence {
    my $self = shift;
    my ($phy) = @_;

    my $alnio = Bio::AlignIO->new(-file => $phy, -format => 'phylip', -longid => 1); 
    my $aln = $alnio->next_aln;
    my $pid = $aln->overall_percentage_identity;
    my $div = 100 - $pid;

    return $div;
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

    perldoc Tephra::LTR::LTRStats


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
