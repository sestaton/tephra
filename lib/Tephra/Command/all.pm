package Tephra::Command::all;
# ABSTRACT: Run all subcommands and generate annotations for all transposon types.

use 5.014;
use strict;
use warnings;
use Pod::Find     qw(pod_where);
use Pod::Usage    qw(pod2usage);
use Cwd           qw(abs_path);
use Capture::Tiny qw(capture_merged);
use File::Spec;
use File::Basename;
use Capture::Tiny;
use Tephra -command;
use Tephra::Config::Reader;
use Tephra::Config::Exe;
use Tephra::Genome::Unmask;
use Tephra::Analysis::Pipeline;
#use Data::Dump::Color;

our $VERSION = '0.13.1';

sub opt_spec {
    return (    
	[ "config|c=s",    "The Tephra configuration file "     ],
	[ "help|h",        "Display the usage menu and exit. "  ],
        [ "man|m",         "Display the full manual. "          ],
    );
}

sub validate_args {
    my ($self, $opt, $args) = @_;

    my $command = __FILE__;
    if ($opt->{man}) {
        system('perldoc', $command) == 0 or die $!;
        exit(0);
    }
    elsif ($opt->{help}) {
        $self->help and exit(0);
    }
    elsif (! $opt->{config}) {
	say STDERR "\n[ERROR]: Required arguments not given.\n";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{config}) {
	say STDERR "\n[ERROR]: The configuration file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $success = _run_all_commands($opt);
}

sub _run_all_commands {
    my ($opt) = @_;

    my (@mask_files, @fas_files, @gff_files, @age_files, @classified_fastas, @unclassified);
    my $config_obj  = Tephra::Config::Reader->new( config => $opt->{config} );
    my $config      = $config_obj->get_configuration;
    my $global_opts = $config_obj->get_all_opts($config);
    #dd $config and exit;

    ## set global options
    my $tephra_obj = Tephra::Analysis::Pipeline->new( global_options => $global_opts, cli_options => $opt );
    my ($tzero, $log) = $tephra_obj->init_tephra($global_opts);

    ## findltrs
    my ($ltr_fas, $ltr_gff) = $tephra_obj->find_ltrs($log);

    ## maskref on LTRs
    my $genome_mask1;
    my $refct = 1;
    my $has_ltrs = 0;
    
    if (defined $ltr_fas && -e $ltr_fas && -s $ltr_fas) {
	$has_ltrs = 1;
	$log->info("Output files - $ltr_gff");
	$log->info("Output files - $ltr_fas");

	# Need to add a 3-letter code ("RLX") to avoid warnings in masking.
	# The LTRs/TRIMs will be classified in a later step and the proper codes will
	# be applied.
	#my $tmp_ltr_ref = $tephra_obj->make_temp_reference_for_masking($ltr_fas, 'LTRs');

	$genome_mask1 = $tephra_obj->mask_reference({
	    log         => $log, 
	    log_results => 0,
	    config      => $config, 
	    reference   => $global_opts->{genome}, 
	    #database    => $tmp_ltr_ref,
	    database    => $ltr_fas,
	    dbtype      => 'LTRs', 
	    ref_count   => $refct });

	if (defined $genome_mask1 && -e $genome_mask1) { 
	    $refct++;
	    $log->info("Output files - $genome_mask1");
	    push @mask_files, $genome_mask1;
	    push @mask_files, $genome_mask1.'.log';
	    #unlink $tmp_ltr_ref;
	}
    }

    ## TRIMs
    my $has_trims = 0;
    my $trims_ref = defined $genome_mask1 && -e $genome_mask1 ? $genome_mask1 : $global_opts->{genome};
    my ($trims_fas, $trims_gff) = $tephra_obj->find_trims($log, $trims_ref);

    if (-e $trims_fas && -e $trims_gff) {
	$has_trims = 1;
	my $unmask = Tephra::Genome::Unmask->new( genome => $global_opts->{genome}, repeatdb => $trims_fas );
	$unmask->unmask_repeatdb;
    }

    ## combine LTRs and TRIMs
    my $ltr_trim_gff;
    if ($has_ltrs && $has_trims) {
	$ltr_trim_gff = $tephra_obj->combine_ltrs_trims($trims_gff, $ltr_gff);
    }
    elsif ($has_ltrs && !$has_trims) {
	$ltr_trim_gff = $ltr_gff;
    }
    elsif (!$has_ltrs && $has_trims) {
	$ltr_trim_gff = $trims_gff;
    }
    else {
	$ltr_trim_gff = undef;
    }

    my ($ltrc_fas, $ltrc_gff, $ltrc_dir);
    if (defined $ltr_trim_gff) {
	($ltrc_fas, $ltrc_gff, $ltrc_dir) = $tephra_obj->classify_ltrs($log, $ltr_trim_gff);

	if (-e $ltrc_gff && -e $ltrc_fas) {
	    $log->info("Output files - $ltrc_gff");
	    $log->info("Output files - $ltrc_fas");
	    push @fas_files, $ltrc_fas;
	    push @gff_files, $ltrc_gff;
	    push @classified_fastas, $ltrc_fas;
	    push @unclassified, $ltr_fas;
	    push @unclassified, $ltr_gff;
	    push @unclassified, $trims_fas;
	    push @unclassified, $trims_gff;
	}
    }

    ## ltrage
    if (defined $ltrc_gff && -e $ltrc_gff && -s $ltrc_gff) {
	my $ltrage_out = $tephra_obj->calculate_ltrage($config, $log, $ltrc_gff, $ltrc_dir);

        if (-e $ltrage_out) {
            $log->info("Output files - $ltrage_out");
            push @age_files, $ltrage_out;
        }
    }

    ## maskref on classified LTRs and TRIMs
    my $genome_mask2;
    if (defined $ltrc_fas && -e $ltrc_fas && -s $ltrc_fas) {
	$genome_mask2 = $tephra_obj->mask_reference({
            log         => $log,
	    log_results => 0,
            config      => $config,
            reference   => $genome_mask1,
            database    => $ltrc_fas,
            dbtype      => 'classified LTRs/TRIMs',
            ref_count   => $refct });

	if (defined $genome_mask2 && -e $genome_mask2) { 
	    $refct++;
	    $log->info("Output files - $genome_mask2");
	    push @mask_files, $genome_mask2;
	    push @mask_files, $genome_mask2.'.log';
	}
    }

    ## sololtrs
    my ($sololtr_gff, $sololtr_rep, $sololtr_fas);
    if (defined $genome_mask1 && -e $genome_mask1 && -e $ltrc_dir) {
	my ($sololtr_gff, $sololtr_rep, $sololtr_fas) = 
	    $tephra_obj->find_sololtrs($config, $log, $genome_mask1, $ltrc_dir);
	
	if (-e $sololtr_gff && -e $sololtr_rep && -e $sololtr_fas) {
	    $log->info("Output files - $sololtr_gff");
	    $log->info("Output files - $sololtr_rep");
	    $log->info("Output files - $sololtr_fas");
	    push @gff_files, $sololtr_gff;
	}
    }

    ## illrecomb
    if (defined $ltrc_fas && -e $ltrc_fas && -s $ltrc_fas) {
	my ($illrec_fas, $illrec_rep, $illrec_stats) = $tephra_obj->find_illrecombination($log, $ltrc_fas);

	if (-e $illrec_fas && -e $illrec_rep && -e $illrec_stats) {
	    $log->info("Output files - $illrec_fas");
	    $log->info("Output files - $illrec_rep");
	    $log->info("Output files - $illrec_stats");
	}
    }

    ## findhelitrons
    my $hel_ref = (defined $genome_mask2 && -s $genome_mask2) ? $genome_mask2 
	        : (defined $genome_mask1 && -s $genome_mask1) ? $genome_mask1 
	        : $global_opts->{genome};

    my ($hel_fas, $hel_gff) = $tephra_obj->find_helitrons($log, $hel_ref);

    if (-e $hel_gff && -e $hel_fas) {
	my $unmask = Tephra::Genome::Unmask->new( genome => $global_opts->{genome}, repeatdb => $hel_fas );
	$unmask->unmask_repeatdb;
	$log->info("Output files - $hel_gff");
	$log->info("Output files - $hel_fas");
	push @fas_files, $hel_fas;
    }

    ## maskref on Helitrons
    my $genome_mask3;
    if (-e $hel_fas && -s $hel_fas) {
	$genome_mask3 = $tephra_obj->mask_reference({
            log         => $log,
	    log_results => 0,
            config      => $config,
            reference   => $hel_ref,
            database    => $hel_fas,
            dbtype      => 'Helitrons',
            ref_count   => $refct });

	if (defined $genome_mask3 && -e $genome_mask3) { 
	    $refct++;
	    $log->info("Output files - $genome_mask3");
	    push @mask_files, $genome_mask3;
	    push @mask_files, $genome_mask3.'.log';
	}
    }

    ## findtirs
    my $tir_ref = (defined $genome_mask3 && -s $genome_mask3) ? $genome_mask3
	        : (defined $genome_mask2 && -s $genome_mask2) ? $genome_mask2 
	        : (defined $genome_mask1 && -s $genome_mask1) ? $genome_mask1 
	        : $global_opts->{genome};

    my ($tir_gff, $tir_fas) = $tephra_obj->find_tirs($log, $tir_ref);

    if (-e $tir_gff && -e $tir_fas) {
	$log->info("Output files - $tir_gff");
	$log->info("Output files - $tir_fas");
    }

    ## classifytirs
    my ($tirc_gff, $tirc_fas, $tirc_dir);
    if (-e $tir_gff && -s $tir_gff) {
	($tirc_fas, $tirc_gff, $tirc_dir) = $tephra_obj->classify_tirs($log, $tir_ref, $tir_gff);

	if (-e $tirc_gff && -e $tirc_fas) {
	    $log->info("Output files - $tirc_gff");
	    $log->info("Output files - $tirc_fas");
	    push @gff_files, $tirc_gff;
	    push @fas_files, $tirc_fas;
	    push @classified_fastas, $tirc_fas;
	    push @unclassified, $tir_fas;
	    push @unclassified, $tir_gff;
	}
    }

    ## tirage
    if (defined $tirc_gff && -e $tirc_gff && -s $tirc_gff) {
	my $tirage_out = $tephra_obj->calculate_tirage($config, $log, $tirc_gff, $tirc_dir);

        if (-e $tirage_out) {
            $log->info("Output files - $tirage_out");
	    push @age_files, $tirage_out;
        }
    }

    ## maskref on TIRs
    my $genome_mask4;
    if (-e $tirc_fas && -s $tirc_fas) {
	$genome_mask4 = $tephra_obj->mask_reference({
            log         => $log,
	    log_results => 0,
            config      => $config,
            reference   => $tir_ref,
            database    => $tirc_fas,
            dbtype      => 'TIRs',
            ref_count   => $refct });

	if (defined $genome_mask4 && -e $genome_mask4) {
	    $refct++;
	    $log->info("Output files - $genome_mask4");
	    push @mask_files, $genome_mask4;
	    push @mask_files, $genome_mask4.'.log';
	}
    }

    ## findnonltrs
    my $nonltr_ref = (defined $genome_mask4 && -s $genome_mask4) ? $genome_mask4
	           : (defined $genome_mask3 && -s $genome_mask3) ? $genome_mask3
	           : (defined $genome_mask2 && -s $genome_mask2) ? $genome_mask2 
	           : (defined $genome_mask1 && -s $genome_mask1) ? $genome_mask1 
	           : $global_opts->{genome};

    my ($nonltr_fas, $nonltr_gff) = $tephra_obj->find_nonltrs($log, $nonltr_ref);

    if (-e $nonltr_gff && -e $nonltr_fas) {
	$log->info("Output files - $nonltr_gff");
	$log->info("Output files - $nonltr_fas");
	push @fas_files, $nonltr_fas;
    }

    ## combine full-length elements for masking and fragment discovery
    my $completeRepDB = $tephra_obj->make_combined_repeatdb({ log => $log, files => \@fas_files, task => 'complete' });
    #my $completeRepDB = $customRepDB =~ s/\.fasta/_complete.fasta/r;
    #move $customRepDB, $completeRepDB or die "\n[ERROR]: move failed: $!\n";

    if (-e $completeRepDB) {
        $log->info("Output files - $completeRepDB");
    }

    ## maskref on customRepDB
    my $final_mask = $tephra_obj->mask_reference({
	  log         => $log,
	  log_results => 1,
	  config      => $config,
	  reference   => $global_opts->{genome},
	  database    => $completeRepDB,
	  dbtype      => 'full-length transposon database',
	  ref_count   => $refct });

    push @mask_files, $final_mask.'.log';

    ## findfragments
    my ($fragments_gff, $fragments_fas) = $tephra_obj->find_fragments($log, $completeRepDB, $final_mask);
    push @fas_files, $fragments_fas;

    if (-e $fragments_gff) {
	$log->info("Output files - $fragments_gff");
	$log->info("Output files - $fragments_fas");
    }
    
    my $customRepDB = $tephra_obj->make_combined_repeatdb({ log => $log, files => \@fas_files, task => 'all' });

    if (-e $customRepDB) {
        $log->info("Output files - $customRepDB");
    }

    ## combine GFF3
    $tephra_obj->combine_gff_files($log, \@gff_files, $hel_gff, $nonltr_gff, $fragments_gff);

    ## calculate family similarity
    $tephra_obj->calculate_global_te_similarity($log, $customRepDB);

    ## combine age files
    my $age_ct = $tephra_obj->combine_age_files($log, \@age_files, \@classified_fastas);

    ## clean up
    my $exe_conf = Tephra::Config::Exe->new->get_config_paths;
    my ($gt, $cleanpp) = @{$exe_conf}{qw(gt cleanpp)};

    $tephra_obj->capture_cmd($cleanpp);
    $tephra_obj->capture_cmd($gt, 'clean');
    my @fais = glob "*.fai";
    unlink @fais;
    unlink @mask_files;
    if (@unclassified) { 
	for my $file (@unclassified) {
	    $log->info("Clean up - Removing unannotated file - $file") if -e $file;
	    unlink $file;
	}
    }

    if ($global_opts->{genome_is_compressed}) {
	unlink $global_opts->{genome};
	$log->info("Clean up - Removing temporary file - $global_opts->{genome}");
    }
    if ($global_opts->{repeatdb_is_compressed}) {
	unlink $global_opts->{repeatdb};
	$log->info("Clean up - Removing temporary file - $global_opts->{repeatdb}");
    }
    if ($global_opts->{genefile_is_compressed}) {
        unlink $global_opts->{genefile};
        $log->info("Clean up - Removing temporary file - $global_opts->{genefile}");
    }

    # Log summary of results
    $tephra_obj->log_interval( $tzero, $log );

    return 1;
}

#
# methods
#
sub help {
    my $desc = capture_merged {
        pod2usage(-verbose => 99, -sections => "NAME|DESCRIPTION", -exitval => "noexit",
		  -input => pod_where({-inc => 1}, __PACKAGE__));
    };
    chomp $desc;
    print STDERR<<END
$desc
USAGE: tephra all [-h] [-m]
    -m --man      :   Get the manual entry for a command.
    -h --help     :   Print the command usage.

Required:
    -c|config     :   The Tephra configuration file.

END
}

1;
__END__

=pod

=head1 NAME
                                                                       
 tephra all - Run all subcommands and generate annotations for all transposon types. 

=head1 SYNOPSIS    

 tephra all -c tephra_pipeline_config.yml

=head1 DESCRIPTION

 This command will run all tephra analysis methods and generate a combined GFF3 of annotations,
 along with a FASTA database and masked genome. A table of results showing global repeat composition 
 will be produced for the input genome.

=head1 AUTHOR 

S. Evan Staton, C<< <evan at evanstaton.com> >>

=head1 REQUIRED ARGUMENTS

=over 2

=item -c, --config

 The Tephra configuration file.

=back

=head1 OPTIONS

=over 2

=item -h, --help

 Print a usage statement. 

=item -m, --man

 Print the full documentation.

=back

=cut
