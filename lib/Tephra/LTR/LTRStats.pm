package Tephra::LTR::LTRStats;

use 5.014;
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
use Cwd             qw(getcwd abs_path);
use Bio::DB::HTS::Kseq;
use Parallel::ForkManager;
use Carp 'croak';
use Try::Tiny;
use Tephra::Alignment::Utils;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::Any',
     'Tephra::Role::Run::PAML';

=head1 NAME

Tephra::LTR::LTRStats - Calculate the age distribution of LTR retrotransposons

=head1 VERSION

Version 0.10.00

=cut

our $VERSION = '0.10.00';
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

has dir => (
    is       => 'ro',
    isa      => 'Path::Class::Dir',
    required => 0,
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
    required => 0,
    coerce   => 1,
);

has clean => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

has all => (
    is       => 'ro',
    isa      => 'Bool',
    required => 0,
    default  => 1,
);

#
# methods
#
sub calculate_ltr_ages {
    my $self = shift;
    my $threads = $self->threads;
    my $outfile = $self->outfile;

    my $args = $self->collect_feature_args;
    
    my $resdir = File::Spec->catdir( abs_path($args->{resdir}), 'divergence_time_stats');
    
    unless ( -d $resdir ) {
	make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $t0 = gettimeofday();
    my $ltrrts  = 0;
    my $logfile = File::Spec->catfile($resdir, 'all_aln_reports.log');
    open my $logfh, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";
    
    my $pm = Parallel::ForkManager->new($threads);
    local $SIG{INT} = sub {
        $log->warn("Caught SIGINT; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      my $t1 = gettimeofday();
			      my $elapsed = $t1 - $t0;
			      my $time = sprintf("%.2f",$elapsed/60);
			      say $logfh basename($ident),
			          " just finished with PID $pid and exit code: $exit_code in $time minutes";
			} );

    for my $type (keys %$args) {
	if ($type eq 'ltrs') {
	    for my $db (@{$args->{$type}{seqs}}) {
		$ltrrts++;
		$pm->start($db) and next;
		$SIG{INT} = sub { $pm->finish };
		$self->process_align_args($db, $resdir);
		
		$pm->finish(0, \$db);
	    }
	}
    }
    $pm->wait_all_children;

    my @agefiles;
    find( sub { push @agefiles, $File::Find::name if -f and /divergence.txt$/ and -s }, $resdir );

    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
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
	find({ wanted => $wanted, preprocess => $process }, $resdir);
	unlink @alnfiles;

	remove_tree( $args->{resdir}, { safe => 1 } );
    }

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $logfh "\n========> Finished calculating insertion time for $ltrrts LTR-RTs in $final_time minutes";
    close $logfh;

    return;
}

sub collect_feature_args {
    my $self = shift;

    my %aln_args;
    if ($self->all || ! $self->dir) {
	my ($files, $wdir) = $self->extract_ltr_features;
        $aln_args{ltrs} = { seqs => $files };
        $aln_args{resdir} = $wdir;
    }
    else {
	my $dir = $self->dir;
	#my $dir = $self->dir->absolute->resolve; 
	my $ltrseqs = $self->_get_exemplar_ltrs($dir);

	if (@$ltrseqs > 0) {
	    $aln_args{ltrs} = { seqs => $ltrseqs };
	    $aln_args{resdir} = $dir;
	}
	else {
	    unless (-e $self->gff) {
		croak "\n[ERROR]: No exemplar files were found in the input directory ".
		    "and the input GFF file does not appear to exist. Exiting.\n";
	    }

	    warn "\n[WARNING]: No exemplar files were found in the input directory. LTR age will be ".
		"calculated from LTR-RT elements in the input GFF.\n";

	    my ($files, $wdir) = $self->extract_ltr_features;
	    $aln_args{ltrs} = { seqs => $files };
	    $aln_args{resdir} = $wdir;
	}
    }

    return \%aln_args;
}

sub extract_ltr_features {
    my $self = shift;
    my $fasta = $self->genome->absolute->resolve;
    my $gff   = $self->gff->absolute->resolve;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $dir = File::Spec->catdir($path, $name.'_ltrages');
    unless ( -d $dir ) {
	make_path( $dir, {verbose => 0, mode => 0771,} );
    }

    my ($header, $features) = $self->collect_gff_features($gff);
    my $index = $self->index_ref($fasta);

    my ($family, %ltrs, %seen, %coord_map);
    for my $rep_region (keys %$features) {
        for my $ltr_feature (@{$features->{$rep_region}}) {
	    if ($ltr_feature->{type} =~ /(?:LTR|TRIM)_retrotransposon/) {
		my $elem_id = @{$ltr_feature->{attributes}{ID}}[0];
		$family  = @{$ltr_feature->{attributes}{family}}[0];
		my ($start, $end) = @{$ltr_feature}{qw(start end)};
		my $key = join "||", $family, $elem_id, $start, $end;
		$ltrs{$key}{'full'} = join "||", @{$ltr_feature}{qw(seq_id type start end)};
		$coord_map{$elem_id} = join "||", @{$ltr_feature}{qw(seq_id start end)};
	    }
	    if ($ltr_feature->{type} eq 'long_terminal_repeat') {
		my $parent = @{$ltr_feature->{attributes}{Parent}}[0];
		my ($seq_id, $pkey) = $self->get_parent_coords($parent, \%coord_map);
		if ($seq_id eq $ltr_feature->{seq_id}) {
		    my ($type, $start, $end, $strand) = 
			@{$ltr_feature}{qw(type start end strand)};
		    $strand //= '?';
		    my $ltrkey = join "||", $seq_id, $type, $start, $end, $strand;
		    my $parent_key = join "||", $family, $pkey;
		    push @{$ltrs{$parent_key}{'ltrs'}}, $ltrkey unless exists $seen{$ltrkey};
		    $seen{$ltrkey} = 1;
		}
	    }
	}
    }

    my @files;
    my $ltrct = 0;
    my $orientation;
    for my $ltr (nsort keys %ltrs) {
	my ($family, $element, $rstart, $rend) = split /\|\|/, $ltr;
	my ($seq_id, $type, $start, $end) = split /\|\|/, $ltrs{$ltr}{'full'};
	my $ltr_file = join "_", $family, $element, $seq_id, $start, $end, 'ltrs.fasta';
	my $ltrs_out = File::Spec->catfile($dir, $ltr_file);
	die "\n[ERROR]: $ltrs_out exists. This will cause problems downstream. Please remove the previous ".
	    "results and try again. Exiting.\n" if -e $ltrs_out;
	push @files, $ltrs_out;
	open my $ltrs_outfh, '>>', $ltrs_out or die "\n[ERROR]: Could not open file: $ltrs_out\n";

	for my $ltr_repeat (@{$ltrs{$ltr}{'ltrs'}}) {
	    #ltr: Contig57_HLAC-254L24||long_terminal_repeat||60101||61950||+
            my ($src, $ltrtag, $s, $e, $strand) = split /\|\|/, $ltr_repeat;
	    my $ltrid;

            if ($ltrct) {
		$orientation = '5prime' if $strand eq '-';
                $orientation = '3prime'  if $strand eq '+';
                $orientation = 'unk-prime-r' if $strand eq '?';
		$ltrid = join "_", $orientation, $family, $element, $src, $s, $e;
		$self->write_element_parts($index, $src, $s, $e, $ltrs_outfh, $ltrid);
                $ltrct = 0;
            }
            else {
		$orientation = '5prime' if $strand eq '+';
                $orientation = '3prime' if $strand eq '-';
                $orientation = 'unk-prime-f' if $strand eq '?';
		$ltrid = join "_", $orientation, $family, $element, $src, $s, $e;
		$self->write_element_parts($index, $src, $s, $e, $ltrs_outfh, $ltrid);
                $ltrct++;
            }
        }
    }

    return (\@files, $dir);
}

sub process_baseml_args {
    my $self = shift;
    my ($phy, $dnd, $resdir) = @_;

    my $cwd = getcwd();
    my ($pname, $ppath, $psuffix) = fileparse($phy, qr/\.[^.]*/);
    my $divfile = File::Spec->catfile( abs_path($ppath), $pname.'-divergence.txt' );
    $divfile = basename($divfile);

    my $utils = Tephra::Alignment::Utils->new;
    my $divergence = $utils->check_divergence($phy);

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
	$element =~ s/_ltrs.*_muscle-out.*//;
	open my $divout, '>', $divfile or die "\n[ERROR]: Could not open divergence file: $divfile\n";
	say $divout join "\t", $element, $divergence , '0', '0';
	close $divout;
	my $dest_file = File::Spec->catfile($resdir, $divfile);
	copy($divfile, $dest_file) or die "\n[ERROR]: Copy failed: $!\n";
	unlink $divfile;
    }
}

sub process_align_args {
    my $self = shift;
    my ($db, $resdir) = @_;

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $pdir = File::Spec->catdir( abs_path($path), $name.'_pamltmp' );
    make_path( $pdir, {verbose => 0, mode => 0771,} );

    my $fas = File::Spec->catfile($pdir, $name.$suffix);
    copy($db, $fas) or die "\n[ERROR]: Copy failed: $!\n";
    unlink $db;

    my $tre  = File::Spec->catfile($pdir, $name.'.dnd');
    my $aln  = File::Spec->catfile($pdir, $name.'_muscle-out.aln');
    my $dnd  = File::Spec->catfile($pdir, $name.'_muscle-out.dnd');
    my $alog = File::Spec->catfile($pdir, $name.'_muscle-out.alnlog');
    my $tlog = File::Spec->catfile($pdir, $name.'_muscle-out.trelog');

    my $muscmd = "muscle -clwstrict -in $fas -out $aln 2>$alog";
    my $trecmd = "muscle -maketree -in $fas -out $tre -cluster neighborjoining 2>$tlog";
    my $status = $self->capture_cmd($muscmd);
    unlink $db && return if $status =~ /failed/i;
    $status = $self->capture_cmd($trecmd);
    unlink $db && return if $status =~ /failed/i;

    my $utils = Tephra::Alignment::Utils->new;
    my $phy = $utils->parse_aln($aln, $tre, $dnd);
    $self->process_baseml_args($phy, $dnd, $resdir);
    
    return;
}

sub _get_exemplar_ltrs {
    my $self = shift;
    my ($dir) = @_;

    my (@dirs, @ltrseqs);
    find( sub { push @dirs, $File::Find::name if -d && /_copia\z|_gypsy\z|_unclassified\z/ }, $dir);
    croak "\n[ERROR]: Could not find the expected sub-directories ending in 'copia', 'gypsy' and 'unclassified'. Please ".
        "check input. Exiting.\n" unless @dirs;

    for my $sfdir (@dirs) {
	my ($ltrfile, %ltrfams);
	find( sub { $ltrfile = $File::Find::name if -f and /exemplar_repeats.fasta$/ }, $sfdir);
	unless (defined $ltrfile) {
	    say STDERR "\n[WARNING]: No exemplar LTR file was found in: $sfdir.";
	    say STDERR "This is likely because there were no families identified by the 'classifyltrs' command for this superfamily.";
	    say STDERR "You can try the 'ltrage' command again with the --all flag to process all LTR-RTs.\n";
	    next;
	}

	my $kseq = Bio::DB::HTS::Kseq->new($ltrfile);
	my $iter = $kseq->iterator();
	
	while ( my $seq = $iter->next_seq() ) {
	    my $id  = $seq->name;
	    my $seq = $seq->seq;
	    if ($id =~ /^[35]prime_(RL[CGX]_family\d+)_(?:LTR|TRIM)_retrotransposon.*/) {
		my $family = $1;
		push @{$ltrfams{$family}}, { id => $id, seq => $seq };
	    }
	}
	
	for my $family (keys %ltrfams) {
	    my $outfile = File::Spec->catfile($dir, $family.'_exemplar_ltrseqs.fasta');
	    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	    for my $pair (@{$ltrfams{$family}}) {
		say $out join "\n", ">".$pair->{id}, $pair->{seq};
	    }
	    close $out;
	    push @ltrseqs, $outfile;
	}
    }

    return \@ltrseqs;
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

    perldoc Tephra::LTR::LTRStats


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
