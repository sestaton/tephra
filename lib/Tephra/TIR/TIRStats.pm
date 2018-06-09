package Tephra::TIR::TIRStats;

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

Tephra::TIR::TIRStats - Calculate the age distribution of TIR transposons

=head1 VERSION

Version 0.11.0

=cut

our $VERSION = '0.11.0';
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
sub calculate_tir_ages {
    my $self = shift;
    my $threads = $self->threads;
    my $outfile = $self->outfile; 

    my $args = $self->collect_feature_args;
    #dd $args; ## debug
    
    my $resdir = File::Spec->catdir($args->{resdir}, 'divergence_time_stats');
    
    unless ( -d $resdir ) {
	make_path( $resdir, {verbose => 0, mode => 0771,} );
    }
    
    my $t0 = gettimeofday();
    my $tirts = 0;
    my $logfile = File::Spec->catfile( abs_path($resdir), 'all_aln_reports.log' );
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
	if ($type eq 'tirs') {
	    for my $db (@{$args->{$type}{seqs}}) {
		$tirts++;
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
    say $out join "\t", "TIR-ID", "Divergence", "Age", "Ts:Tv";

    for my $file (@agefiles) {
	$self->collate($file, $out);
    }
    close $out;

    if ($self->clean) {
	my @alnfiles;
	my $wanted  = sub { push @alnfiles, $File::Find::name 
				if -f && /tirs.fasta$|\.aln$|\.log$|\.phy$|\.dnd$|\.ctl$|\-divergence.txt$/ };
	my $process = sub { grep ! -d, @_ };
	find({ wanted => $wanted, preprocess => $process }, $resdir);
	unlink @alnfiles;

	remove_tree( $args->{resdir}, { safe => 1 } );
    }

    my $t2 = gettimeofday();
    my $total_elapsed = $t2 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);

    say $logfh "\n========> Finished calculating insertion time for $tirts TIRs in $final_time minutes";
    close $logfh;

    return;
}

sub collect_feature_args {
    my $self = shift;

    my (@tirs, %aln_args);
    if ($self->all || ! $self->dir) {
	my ($files, $wdir) = $self->extract_tir_sequences;
        $aln_args{tirs} = { seqs => $files };
        $aln_args{resdir} = $wdir;
    }
    else {
	my $dir = $self->dir->absolute->resolve;
	my $wanted  = sub { push @tirs, $File::Find::name if -f && /exemplar_repeats.fasta$/ };
	my $process = sub { grep ! -d, @_ };
	find({ wanted => $wanted, preprocess => $process }, $dir);

	if (@tirs > 0) {
	    $aln_args{tirs} = { seqs => \@tirs };
	    $aln_args{resdir} = $dir;
	}
	else {
	    unless (-e $self->gff) {
		croak "\n[ERROR]: No exemplar files were found in the input directory ".
		    "and the input GFF file does not appear to exist. Exiting.\n\n";
	    }

	    warn "\n[WARNING]: No exemplar files were found in the input directory. TIR age will be ".
		"calculated from TIR elements in the input GFF.\n\n";

	    my ($files, $wdir) = $self->extract_tir_sequences;
	    $aln_args{tirs} = { seqs => $files };
	    $aln_args{resdir} = $wdir;
	}
    }

    return \%aln_args;
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
	$element =~ s/_tirs_muscle-out.*//;
	open my $divout, '>', $divfile or die "\n[ERROR]: Could not open divergence file: $divfile\n";
	say $divout join "\t", $element, $divergence , '0', '0';
	close $divout;
	my $dest_file = File::Spec->catfile( abs_path($resdir), $divfile );
	copy($divfile, $dest_file) or die "\n[ERROR]: Copy failed: $!";
	unlink $divfile;
    }
}

sub process_align_args {
    my $self = shift;
    my ($db, $resdir) = @_;

    my $seqct = $self->_check_tirct($db);
    unlink $db && return unless $seqct == 2;

    my ($name, $path, $suffix) = fileparse($db, qr/\.[^.]*/);
    my $pdir = File::Spec->catdir( abs_path($path), $name.'_pamltmp' );
    make_path( $pdir, {verbose => 0, mode => 0771,} );
	
    my $fas = File::Spec->catfile( abs_path($pdir), $name.$suffix );
    copy($db, $fas) or die "\n[ERROR]: Copy failed: $!";
    my $tre  = File::Spec->catfile( abs_path($pdir), $name.'.dnd' );
    my $aln  = File::Spec->catfile( abs_path($pdir), $name.'_muscle-out.aln' );
    my $dnd  = File::Spec->catfile( abs_path($pdir), $name.'_muscle-out.dnd' );
    my $alog = File::Spec->catfile( abs_path($pdir), $name.'_muscle-out.alnlog' );
    my $tlog = File::Spec->catfile( abs_path($pdir), $name.'_muscle-out.trelog' );
    
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

sub write_tir_parts {
    my $self = shift;
    my ($index, $loc, $elem, $start, $end, $out, $orient, $family) = @_;

    my ($seq, $length) = $self->get_full_seq($index, $loc, $start, $end);

    # need to reverse-complement the inverted seq
    my $utils = Tephra::Alignment::Utils->new;
    $seq = $utils->revcom($seq) if $orient =~ /3prime|prime-r/;

    my $id;
    $id = join "_", $family, $elem, $loc, $start, $end if !$orient;
    $id = join "_", $orient, $family, $elem, $loc, $start, $end if $orient; # for unique IDs with clustalw

    #$self->write_element_parts($index, $loc, $start, $end, $out, $id);
    say $out join "\n", ">".$id, $seq;

    return;
}

sub _check_tirct {
    my $self = shift;
    my ($db) = @_;
    
    my $kseq = Bio::DB::HTS::Kseq->new($db);
    my $iter = $kseq->iterator;

    my $ct = 0;
    while (my $seq = $iter->next_seq) {
	$ct++ if defined $seq->seq;
    }

    return $ct;
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

    perldoc Tephra::TIR::TIRStats


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
