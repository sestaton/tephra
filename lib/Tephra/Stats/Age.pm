package Tephra::Stats::Age;

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
use Tephra::Config::Exe;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::File',
     'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::Any',
     'Tephra::Role::Run::PAML',
     'Tephra::LTR::Role::Utils',
     'Tephra::TIR::Role::Utils';

=head1 NAME

Tephra::Stats::Age - Calculate the age distribution of LTR/TIR transposons

=head1 VERSION

Version 0.13.0

=cut

our $VERSION = '0.13.0';
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

has type => (
      is        => 'ro',
      isa       => 'Maybe[Str]',
      predicate => 'has_type',
      required  => 0,
);
#
# methods
#
sub calculate_ages {
    my $self = shift;
    my $threads = $self->threads;
    my $outfile = $self->outfile;

    my $args = $self->collect_feature_args;
    #dd $args and exit;
    my $resdir = File::Spec->catdir( abs_path($args->{resdir}), 'divergence_time_stats');
    
    unless ( -d $resdir ) {
	make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $t0 = gettimeofday();
    my $tes  = 0;
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
	if ($type eq 'repeats') {
	    for my $db (nsort @{$args->{$type}{seqs}}) {
		$tes++;
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
    say $out join "\t", "ID", "Divergence", "Age", "Ts:Tv";

    for my $file (nsort @agefiles) {
	$self->collate($file, $out);
    }
    close $out;

    if ($self->clean) {
	my @alnfiles;
	my $wanted  = sub { push @alnfiles, $File::Find::name 
				if -f && /(?:ltrs|tirs).fasta$|\.aln$|\.log$|\.phy$|\.dnd$|\.ctl$|\-divergence.txt$/ };
	my $process = sub { grep ! -d, @_ };
	find({ wanted => $wanted, preprocess => $process }, $resdir);
	unlink @alnfiles;

	remove_tree( $args->{resdir}, { safe => 1 } );
    }

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $order = $self->type =~ /ltr/i ? 'LTR-RTs' : 'TIRs';

    say $logfh "\n========> Finished calculating insertion time for $tes $order transposons in $final_time minutes";
    close $logfh;

    return;
}

sub collect_feature_args {
    my $self = shift;
    my $type = uc($self->type);
    my $gff  = $self->gff;

    ## The control statement below checks if we want to:
    ## 1) calculate age on all elements, or
    ## 2) only calculate age on exemplars
    ## NB: the third part of this structure would handle the case where no exemplars are
    ## found and the '--all' flag was not set. In that case, age on all elements would
    ## be calculated instead of exiting. Currently (v0.11.1), the loop would never get
    ## to that point because the GFF3 file and the '--all' flag, or the input directory
    ## must given or the program with halt. It would be easy to allow this behavior in
    ## in a future release once these methods are deemed stable. At the time of a major 
    ## refactoring of these methods, it makes more sense to limit the options.
    my (%aln_args, @seqs);
    if ($self->all && $self->gff) {
	my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
	my $wdir = $type =~ /ltr/i ? File::Spec->catdir($path, $name.'_ltrages') : File::Spec->catdir($path, $name.'_tirages');
	unless ( -d $wdir ) {
	    make_path( $wdir, {verbose => 0, mode => 0771,} );
	}
	
	my ($files, $swdir) = $type =~ /ltr/i ? $self->extract_ltr_sequences : $self->extract_tir_sequences;
	$aln_args{repeats} = { seqs => $files };
	$aln_args{resdir} = $swdir;
    }
    elsif (! $self->all && $self->dir) { 
	my $dir = $self->dir;
	my ($name, $path, $suffix) = fileparse($dir, qr/\.[^.]*/);

	my $wdir = $type =~ /ltr/i ? File::Spec->catdir($path, $name.'_ltrages') : File::Spec->catdir($path, $name.'_tirages');
	unless ( -d $wdir ) {
	    make_path( $wdir, {verbose => 0, mode => 0771,} );
	}

	if ($type =~ /ltr/i) { 
	    my $ltrseqs = $self->get_exemplar_ltrs_for_age($dir, $wdir);
	    @seqs = @$ltrseqs;
	}
	elsif ($type =~ /tir/i) {
	    my $tirseqs = $self->get_exemplar_tirs_for_age($dir, $wdir);
	    @seqs = @$tirseqs;
	}

	if (@seqs > 0) {
	    $aln_args{repeats} = { seqs => \@seqs };
	    $aln_args{resdir} = $wdir; 
	}
    }
    elsif ($self->gff) {
	unless (-e $self->gff) {
	    croak "\n[ERROR]: No exemplar files were found in the input directory ".
		"and the input GFF file does not appear to exist. Exiting.\n";
	}
	
	warn "\n[WARNING]: No exemplar files were found in the input directory. $type age will be ".
	    "calculated from $type elements in the input GFF.\n";
	
	my ($files, $wdir) = $type =~ /ltr/i ? $self->extract_ltr_sequences : $self->extract_tir_sequences;
	$aln_args{repeats} = { seqs => $files };
	$aln_args{resdir} = $wdir;
    }

    #dd \%aln_args and exit;
    return \%aln_args;
}

sub process_baseml_args {
    my $self = shift;
    my ($phy, $dnd, $resdir) = @_;

    my $cwd = getcwd();
    my ($pname, $ppath, $psuffix) = fileparse($phy, qr/\.[^.]*/);
    my $divfile = File::Spec->catfile( abs_path($ppath), $pname.'-divergence.txt' );
    my $divfile_name = basename($divfile);

    my $utils = Tephra::Alignment::Utils->new;
    my $divergence = $utils->check_divergence($phy);

    if ($divergence > 0) {
	my $baseml_args = $self->create_baseml_files({ phylip => $phy, treefile => $dnd });
	$self->run_baseml({ working_dir => $ppath, current_dir => $cwd }); 
	$self->parse_baseml({ divergence_file => $divfile_name,
			      results_dir     => $resdir,
			      treefile        => $dnd,
			      phylip          => $phy,
			      outfile         => $baseml_args->{outfile},
			      control_file    => $baseml_args->{control_file} });
    }
    else {
	my $element = basename($phy);
	$element =~ s/_(?:ltrs|tirs).*_muscle-out.*//;
	open my $divout, '>', $divfile or die "\n[ERROR]: Could not open divergence file: $divfile\n";
	say $divout join "\t", $element, $divergence , '0', '0';
	close $divout;
	my $dest_file = File::Spec->catfile($resdir, $divfile_name);

	copy($divfile, $dest_file) or die "\n[ERROR]: Copy failed for $divfile to $dest_file WD:$cwd: $!\n";
	unlink $divfile;
	remove_tree( $ppath, { safe => 1 } );
    }

    return;
}

sub process_align_args {
    my $self = shift;
    my ($db, $resdir) = @_;

    my $config  = Tephra::Config::Exe->new->get_config_paths;
    my ($muscle) = @{$config}{qw(muscle)};

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $pdir = File::Spec->catdir( abs_path($path), $name.'_pamltmp' );
    make_path( $pdir, {verbose => 0, mode => 0771,} );

    my $fas = File::Spec->catfile($pdir, $name.$suffix);
    copy($db, $fas) or die "\n[ERROR]: Copy failed before alignment for '$path.$name.$suffix': $!\n";
    unlink $db;

    my $tre  = File::Spec->catfile($pdir, $name.'.dnd');
    my $aln  = File::Spec->catfile($pdir, $name.'_muscle-out.aln');
    my $dnd  = File::Spec->catfile($pdir, $name.'_muscle-out.dnd');
    my $alog = File::Spec->catfile($pdir, $name.'_muscle-out.alnlog');
    my $tlog = File::Spec->catfile($pdir, $name.'_muscle-out.trelog');

    my $muscmd = "$muscle -clwstrict -in $fas -out $aln 2>$alog";
    my $trecmd = "$muscle -maketree -in $fas -out $tre -cluster neighborjoining 2>$tlog";
    my $status = $self->capture_cmd($muscmd);
    unlink $db && return if $status =~ /failed/i;
    $status = $self->capture_cmd($trecmd);
    unlink $db && return if $status =~ /failed/i;

    my $utils = Tephra::Alignment::Utils->new;
    my $phy = $utils->parse_aln($aln, $tre, $dnd);
    if (-s $phy) { 
	$self->process_baseml_args($phy, $dnd, $resdir);
    }
    else {
	say STDERR "\n[ERROR]: Phylip conversion of '$phy' did not complete. Check '$aln' and '$tre' for debugging. Exiting.\n";
	exit(1);
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

    perldoc Tephra::Stats::Age


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
