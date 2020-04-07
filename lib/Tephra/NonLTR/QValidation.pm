package Tephra::NonLTR::QValidation;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Parallel::ForkManager;
use IPC::System::Simple qw(system capture);
use File::Path          qw(make_path);
use Cwd                 qw(abs_path);
use Time::HiRes         qw(gettimeofday);
use Bio::SearchIO;
use File::Find;
use File::Spec;
use File::Basename;
use Carp 'croak';
#use Data::Dump::Color;
use Tephra::Config::Exe;
use namespace::autoclean;

with 'Tephra::NonLTR::Role::PathUtils';

=head1 NAME

Tephra::NonLTR::QValidation - Valididate non-LTR search results against all models (adapted from MGEScan-nonLTR)

=head1 VERSION

Version 0.12.6

=cut

our $VERSION = '0.12.6';
$VERSION = eval $VERSION;

has fasta   => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has outdir  => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has phmmdir => ( is => 'ro', isa => 'Path::Class::Dir',  required => 1, coerce => 1 );
has threads => ( is => 'ro', isa => 'Int',  predicate => 'has_threads', lazy => 1, default => 1 );

sub validate_q_score {
    my $self = shift;
    my $dir     = $self->outdir->absolute->resolve;
    my $hmm_dir = $self->phmmdir->absolute->resolve;
    my $genome  = $self->fasta->absolute->resolve;

    my $t0 = gettimeofday();
    my $hmmsearch = $self->find_hmmsearch;
    my @en_clade  = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1');
    my @all_clade = (@en_clade, 'R2', 'CRE');
    my $tree_dir;

    my $infodir = File::Spec->catdir($dir, 'info');
    my $fulldir = File::Spec->catdir($dir, 'info', 'full');
    make_path( $infodir, {verbose => 0, mode => 0771,} );
    make_path( $fulldir, {verbose => 0, mode => 0771,} );
    my $logfile = File::Spec->catfile($fulldir, 'tephra_non-ltrs_qvalidation.log');
    open my $log, '>>', $logfile or die "\n[ERROR]: Could not open file: $logfile\n";

    $self->get_full_frag($genome, $dir, \@all_clade);
    
    # get Q value after running pHMM for EN in full elements
    my $validation_dir = File::Spec->catdir($dir, 'info', 'validation');
    make_path( $validation_dir, {verbose => 0, mode => 0771,} );
    
    my @cladedirs = map { File::Spec->catdir($fulldir, $_) } @all_clade;

    # get domain seq
    my $pm = Parallel::ForkManager->new(2);
    local @{$SIG}{qw(INT TERM)} = sub {
        $log->warn("Caught SIGINT or SIGTERM; Waiting for child processes to finish.");
        $pm->wait_all_children;
        exit 1;
    };

    $pm->run_on_finish( sub { my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_ref) = @_;
			      my $domain = uc($data_ref->{domain});

			      my $t1 = gettimeofday();
                              my $elapsed = $t1 - $t0;
                              my $time = sprintf("%.2f",$elapsed/60);
                              say $log "QValidataion for $domain just finished with PID $pid and exit code: ".
				  "$exit_code in $time minutes";
                        } );

    for my $domain (qw(en rt)) { 
	$pm->start($domain) and next;
	@{$SIG}{qw(INT TERM)} = sub { $pm->finish };
	my $dom = $self->get_domain_coords({ domain     => $domain,
					     dir        => $dir,
					     hmm_dir    => $hmm_dir,
					     en_clade   => \@en_clade,
					     all_clade  => \@all_clade,
					     clade_dirs => \@cladedirs,
					     valid_dir  => $validation_dir,
					     genome     => $genome });

	#unless (defined $dom->{domain}) {
	    #say STDERR "\n[MAYDAY]: domain => $domain, $dir, $hmm_dir, $validation_dir, $genome";
	#}

	$pm->finish(0, $dom);
    }
    $pm->wait_all_children;
    close $log;

    return;
}

sub get_domain_coords {
    my $self = shift;
    my ($obj) = @_;

    my $clade = $obj->{domain} eq 'en' ? $obj->{en_clade} : $obj->{all_clade};

    $self->get_domain_for_full_frag($obj, $clade);
    $self->get_QVal($obj, $clade);

    return { domain => $obj->{domain} };
}

sub get_QVal {
    my $self = shift;
    my ($obj, $dom_clade) = @_;

    my ($domain, $clade_dirs, $validation_dir, $hmm_dir) =
        @{$obj}{qw(domain clade_dirs valid_dir hmm_dir)};

    my $validation_file = File::Spec->catfile($validation_dir, $domain);
    my $evalue_file     = File::Spec->catfile($validation_dir, $domain.'_evalue');

    for my $clade (@$clade_dirs) {
        my $name = basename($clade);
        my $seq  = File::Spec->catfile($clade, $name.'.'.$domain.'.pep');
        if (-e $seq) {
            $self->vote_hmmsearch($seq, $hmm_dir, $domain, $validation_file, $evalue_file, $dom_clade);
        }
    }

    return;
}

sub get_domain_for_full_frag {
    my $self = shift;
    my ($obj, $clade_dir) = @_;

    my ($genome, $domain, $dir, $hmm_dir) =
	@{$obj}{qw(genome domain dir hmm_dir)};

    for my $clade (@$clade_dir) {
	my $resdir   = File::Spec->catdir($dir, 'info', 'full', $clade);
	my $pep_file = File::Spec->catfile($resdir, $clade.'.pep');
	my $dna_file = File::Spec->catfile($resdir, $clade.'.dna');
       
	if (-e $pep_file) {
	    my $phmm_file       = File::Spec->catfile($hmm_dir, $clade.'.'.$domain.'.hmm');
	    my $result_pep_file = File::Spec->catfile($resdir, $clade.'.'.$domain.'.pe');
	    my $result_dna_file = File::Spec->catfile($resdir, $clade.'.'.$domain.'.dna');
	    
	    my $flag = 2;  #1: protein-protein, 2: protein-dna 
	    $self->get_domain_pep_seq($pep_file, $phmm_file, $result_pep_file);
	    $self->get_domain_dna_seq($pep_file, $phmm_file, $result_dna_file, $dna_file, $flag);
	    $self->_add_clade_to_header($clade, $result_pep_file);
	}
    }

    return;
}

sub get_full_frag {
    my $self = shift;
    my ($genome, $dir, $all_clade) = @_;

    for my $clade (@$all_clade) {
	# create a clade dir
	my $clade_dir = File::Spec->catdir($dir, 'info', 'full', $clade);
	my $file_f    = File::Spec->catdir($dir, 'f', 'out2', $clade.'_full');
	my $file_b    = File::Spec->catdir($dir, 'b', 'out2', $clade.'_full');

	if (-e $file_f || -e $file_b){
	    make_path( $clade_dir, {verbose => 0, mode => 0771,} );
	}

	# copy full length in + strand
	if (-e $file_f) {
	    my $of = File::Spec->catfile($clade_dir, $clade.'.dna');
	    open my $fh_out, '>', $of or die"\n[ERROR]: Could not open file: $of";
	    $self->_add_strand_to_header($file_f, $fh_out, '+');
	    close $fh_out;
	}

	# copy full length in - strand
	if (-e $file_b) {
	    my $of = File::Spec->catfile($clade_dir, $clade.'.dna');
            open my $fh_out, '>>', $of or die"\n[ERROR]: Could not open file: $of";
            $self->_add_strand_to_header($file_b, $fh_out, '-');
	    close $fh_out;
	}

	# translate
	my $dna_file = File::Spec->catfile($clade_dir, $clade.'.dna');
	my $pep_file = File::Spec->catfile($clade_dir, $clade.'.pep');   
	if (-e $dna_file){
	    my $config = Tephra::Config::Exe->new->get_config_paths;
	    my ($transeq) = @{$config}{qw(transeq)};
	    system([0..5], $transeq, '-frame', 'F', '-sequence', $dna_file, '-outseq', $pep_file, '-trim', 'yes', '-auto');
	}
    }

    return;
}

sub get_domain_pep_seq {
    my $self = shift;
    my ($pep_file, $phmm_file, $result_pep_file) = @_;

    my $hmmsearch = $self->find_hmmsearch;
    my @hmm_results = capture([0..5], $hmmsearch, $phmm_file, $pep_file);
    my (%domain_start, %domain_end, %result_start, %result_end, %uniq_head);

    my ($name, $path, $suffix) = fileparse($pep_file, qr/\.[^.]*/);
    my $hmmout = File::Spec->catfile($path, $name.'_hmmsearch.txt');
    open my $o, '>', $hmmout or die "\n[ERROR]: Could not open file: $hmmout";
    print $o @hmm_results;
    close $o;
    my $hmmer_in = Bio::SearchIO->new(-file => $hmmout, -format => 'hmmer');

    while ( my $result = $hmmer_in->next_result ) {    
        while ( my $hit = $result->next_hit ) {
	    my $hitid = $hit->name;
            while ( my $hsp = $hit->next_hsp ) {
		my $hstart = $hsp->start('hit');
		my $hend   = $hsp->end('hit');
		#say "DEBUG: $hmmout => hmmermatch qvalidation: $1";
		#Ha1.fa_79699679_3    1/1     668   899 ..     1   260 []   185.0  4.4e-55
		#my @temp = split /\s+/, $1;
		my $key      = substr($hitid, 0, length($hitid));
		my $uniq_key = substr($hitid, 0, length($hitid)-2);
		#say STDERR join q{ }, "key: $key", "uniq_key: $uniq_key";

		if (not exists $uniq_head{$uniq_key}){
		    $uniq_head{$uniq_key} = 1;
		    $result_start{$key}   = $hstart;
		    $result_end{$key}     = $hend;
		}
	    }
	}
    }

    my $flag = 0;
    my ($head, $seq);

    open my $in, '<', $pep_file or die "\n[ERROR]: Could not open file: $pep_file";
    open my $out, '>', $result_pep_file or die "\n[ERROR]: Could not open file: $result_pep_file";

    while (my $line = <$in>){
	chomp $line;
	if ($line =~ /\>/){
	    if (defined $head && length($head) > 0 && $flag == 1) {
		say $out '>'.$head.'_'.$result_start{$head}.'_'.$result_end{$head}; #put region in header
		say $out substr($seq, $result_start{$head}, eval($result_end{$head}-$result_start{$head}+1));
	    }

	    my @temp = split /\s+/, $line;

	    if (exists $result_start{substr($temp[0], 1, length($temp[0])-1)}) {
		$flag = 1;
		$head = substr($temp[0], 1, length($temp[0])-1);
	    }
	    else {
		$flag = 0;
	    }
	    $seq = "";
	}
	else {
	    if ($flag == 1) {
	        $seq .= $line;
	    }
	}
    }

    if ($flag == 1) {
	say $out '>'.$head.'_'.$result_start{$head}.'_'.$result_end{$head}; # get region in header
	say $out substr($seq, $result_start{$head}, eval($result_end{$head}-$result_start{$head}+1));
    }
    close $in;
    close $out;

    return;
}


sub get_domain_dna_seq {
    my $self = shift;
    my ($pep_file, $phmm_file, $result_dna_file, $dna_file, $flag) = @_;

    my $threads = $self->threads;
    my $hmmsearch = $self->find_hmmsearch;
    my @hmm_results = capture([0..5], $hmmsearch, '--cpu', $threads, $phmm_file, $pep_file);
    my (%domain_start, %domain_end, %result_start, %result_end, %uniq_head);

    #say "debug2 hmmermatch qvalidation: $1";
    #Ha1.fa_79699679_3    1/1     668   899 ..     1   260 []   185.0  4.4e-55
    #my @temp = split /\s+/, $1;
    my ($name, $path, $suffix) = fileparse($pep_file, qr/\.[^.]*/);
    my $hmmout = File::Spec->catfile($path, $name.'_hmmsearch.txt');
    open my $o, '>', $hmmout or die "\n[ERROR]: Could not open file: $hmmout";;
    print $o @hmm_results;
    close $o;
    my $hmmer_in = Bio::SearchIO->new(-file => $hmmout, -format => 'hmmer');

    while ( my $result = $hmmer_in->next_result ) {    
        while ( my $hit = $result->next_hit ) {
            my $hitid = $hit->name;
            while ( my $hsp = $hit->next_hsp ) {
                my $hstart   = $hsp->start('hit');
                my $hend     = $hsp->end('hit');
		my $key      = substr($hitid, 0, length($hitid));
		my $uniq_key = substr($hitid, 0, length($hitid)-2);

		if (not exists $result_start{$uniq_key}) {
		    $uniq_head{$uniq_key} = 1;
		    if ($flag == 1) {
			$result_start{$key} = $hstart;
			$result_end{$key}   = $hend;
		    }
		    elsif ($flag == 2) {
			$result_start{$uniq_key} = $hstart;
			$result_end{$uniq_key}   = $hend;
		    }
		}
	    }
	}
    }

    $flag = 0;
    my $head;
    my $seq;

    open my $in, '<', $dna_file or die "\n[ERROR]: Could not open file: $dna_file";
    open my $out, '>', $result_dna_file or die "\n[ERROR]: Could not open file: $result_dna_file";

    while (my $each_line = <$in>) {
	chomp $each_line;
	if ($each_line =~ /\>/) {
	    if (defined $head && length($head) > 0 && $flag == 1){
		say $out '>'.$head.'_'.$result_start{$head}.'_'.$result_end{$head}; # header
		say $out substr($seq, $result_start{$head}*3-3, 
				eval(($result_end{$head}-$result_start{$head}+1)*3+3));
	    }

	    my @temp = split /\s+/, $each_line;

	    if (exists $result_start{substr($temp[0], 1, length($temp[0])-1)}) {
		$flag = 1;
		$head = substr($temp[0], 1, length($temp[0])-1);
	    }
	    else{
		$flag = 0;
	    }
	    $seq = "";
	}
	else {
	    if ($flag == 1) {
	        $seq .= $each_line;
	    }
	}
    }

    if ($flag == 1) {
	say $out '>'.$head.'_'.$result_start{$head}.'_'.$result_end{$head};
	say $out substr($seq, $result_start{$head}*3-3, eval(($result_end{$head}-$result_start{$head}+1)*3+3));
    }
    close $in;
    close $out;

    return;
}

sub vote_hmmsearch {
    my $self = shift;
    my ($seq, $hmm_dir, $domain, $validation_file, $evalue_file, $en_clade) = @_;

    my $hmmsearch = $self->find_hmmsearch;
    my %evalue;
    my %save_evalue;
    my %clade;
    my %sig;
    my $i;
    my $anno_clade; 

    open my $in, '<', $seq or die "\n[ERROR]: Could not open file: $seq";
    while (my $line = <$in>) {
	if ($line =~ /\>/) {
	    chomp $line;
	    my $uniq_key = substr($line, 1, length($line)-1);
	    $evalue{$uniq_key}      = 1000;
	    $save_evalue{$uniq_key} = 1000;
	    $clade{$uniq_key}       = '-';
	    $sig{$uniq_key}         = 1;
	}
    }
    close $in;
    
    for my $clade (@$en_clade) {
	my $phmm_file = File::Spec->catfile($hmm_dir, $clade.'.'.$domain.'.hmm');
	my @hmm_results = capture([0..5], $hmmsearch, $phmm_file, $seq);

	my ($name, $path, $suffix) = fileparse($seq, qr/\.[^.]*/);
	my $hmmout = File::Spec->catfile($path, $name.'_hmmsearch.txt');

	open my $o, '>', $hmmout or die "\n[ERROR]: Could not open file: $hmmout";;
	print $o @hmm_results;
	close $o;

	my $hmmer_in = Bio::SearchIO->new(-file => $hmmout, -format => 'hmmer');

	while ( my $result = $hmmer_in->next_result ) {    
	    while ( my $hit = $result->next_hit ) {
		my $hitid = $hit->name;
		while ( my $hsp = $hit->next_hsp ) {
		    my $hstart = $hsp->start('hit');
		    my $hend   = $hsp->end('hit');
		    my $e_val  = $hsp->evalue;
		    
		    #say "debug3 hmmermatch qvalidation: $1";
		    #Ha1.fa_79699679_3    1/1     668   899 ..     1   260 []   185.0  4.4e-55
		    #my @temp = split /\s+/, $1;
		    my $uniq_key = substr($hitid, 0, length($hitid));	    
		    #say STDERR "vote_hmmsearch: $uniq_key";
		    next unless (defined $uniq_key && defined $e_val && defined $clade && defined $save_evalue{$uniq_key});
		    # instead of a regex (as in mgescan), we check if variables are defined

		    $save_evalue{$uniq_key} = join "\t", $save_evalue{$uniq_key}, $clade, $e_val;
		    
		    if ($evalue{$uniq_key} > $e_val) {
			$sig{$uniq_key}    = $e_val/$evalue{$uniq_key};
			$evalue{$uniq_key} = $e_val;
			$clade{$uniq_key}  = $clade;
		    }
		    elsif ($evalue{$uniq_key}/$e_val > $sig{$uniq_key}) {
			$sig{$uniq_key} = $evalue{$uniq_key}/$e_val;
		    }
		}
	    }
	}
    }

    open my $out, '>>', $validation_file or die "\n[ERROR]: Could not open file: $validation_file";
    if ($seq =~ /\/((\w|\d)+)\./) {
	$anno_clade = $1;
    }
    say $out "$anno_clade-------------------------------------";
    for my $key (keys %evalue) {
        print $out join "\t", $key, $clade{$key}, $evalue{$key}."\t";
	printf $out "%.1e\n", $sig{$key};
    }
    close $out;

    return;
}

sub _add_strand_to_header {
    my $self = shift;
    my ($file_in, $fh_out, $strand) = @_;
    open my $fh_in, '<', $file_in or die "\n[ERROR]: Could not open file: $file_in\n";

    while (my $line = <$fh_in>) {
	chomp $line;
	if ($line =~ /^>/) {
	    $line .= "_$strand";
	    say $fh_out $line;
	}
	else {
	    say $fh_out $line;
	}
    }
    close $fh_in;

    return;
}

sub _add_clade_to_header {
    my $self = shift;
    my ($clade, $result_file) = @_;

    my $outfile = $result_file.'p';
    open my $in, '<', $result_file or die "\n[ERROR]: Could open file: $result_file\n";
    open my $out, '>', $outfile or die "\n[ERROR]: Could open file: $outfile\n";

    while (my $line = <$in>) {
	chomp $line;
	if ($line =~ /^>/) {
	    $line =~ s/>/>${clade}_/;
	    say $out $line;
	}
	else {
	    say $out $line;
	}
    }
    close $in;
    close $out;

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

    perldoc Tephra::NonLTR::QValidation


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
