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
use Tephra::Analysis::Pipeline;
#use Data::Dump::Color;

our $VERSION = '0.09.3';

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
	say STDERR "\nERROR: Required arguments not given.\n";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{config}) {
	say STDERR "\nERROR: The configuration file does not exist. Check arguments.\n";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $gff = _run_all_commands($opt);
}

sub _run_all_commands {
    my ($opt) = @_;

    my (@mask_files, @fas_files, @gff_files, @age_files, @classified_fastas);
    my $config_obj  = Tephra::Config::Reader->new( config => $opt->{config} );
    my $config      = $config_obj->get_configuration;
    my $global_opts = $config_obj->get_all_opts($config);
    #dd $config and exit;

    ## set global options
    my $tephra_obj = Tephra::Analysis::Pipeline->new( global_options => $global_opts, cli_options => $opt );
    my ($tzero, $log) = $tephra_obj->init_tephra($global_opts);
    #my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);   

    ## findltrs
    my ($ltr_fas, $ltr_gff) = $tephra_obj->find_ltrs($log);

    ## maskref on LTRs
    my $genome_mask1;
    my $refct = 1;
    my $has_ltrs = 0;
    
    if (defined $ltr_fas && -e $ltr_fas && -s $ltr_fas) {
	$has_ltrs = 1;
	$log->info("Output files - $ltr_gff");

	# Need to add a 3-letter code ("RLX") to avoid warnings in masking.
	# The LTRs/TRIMs will be classified in a later step and the proper codes will
	# be applied.
	my $tmp_ltr_ref = $tephra_obj->make_temp_reference_for_masking($ltr_fas, 'LTRs');

	$genome_mask1 = $tephra_obj->mask_reference({
	    log            => $log, 
	    config         => $config, 
	    reference      => $global_opts->{genome}, 
	    database       => $tmp_ltr_ref, 
	    dbtype         => 'LTRs', 
	    ref_count      => $refct });

	if (defined $genome_mask1 && -e $genome_mask1) { 
	    $refct++;
	    $log->info("Output files - $genome_mask1");
	    push @mask_files, $genome_mask1;
	    push @mask_files, $genome_mask1.'.log';
	    unlink $tmp_ltr_ref;
	}
    }

    ## TRIMs
    my $has_trims = 0;
    my $trims_ref = defined $genome_mask1 && -e $genome_mask1 ? $genome_mask1 : $global_opts->{genome};
    my ($trims_fas, $trims_gff) = $tephra_obj->find_trims($log, $trims_ref);

    if (-e $trims_fas && -e $trims_gff) {
	$has_trims = 1;
    }

    ## combine LTRs and TRIMs
    my ($ltrc_gff, $ltrc_fas, $ltrc_dir, $ltr_trim_gff);
    if ($has_ltrs && $has_trims) {
	$ltr_trim_gff = $tephra_obj->combine_ltrs_trims($trims_gff, $ltr_gff);

	## classifyltrs on combined LTRs and TRIMs
	($ltrc_fas, $ltrc_gff, $ltrc_dir) = $tephra_obj->classify_ltrs($log, $ltr_trim_gff);

	if (-e $ltrc_gff && -e $ltrc_fas) {
	    $log->info("Output files - $ltrc_gff");
	    $log->info("Output files - $ltrc_fas");
	    push @fas_files, $ltrc_fas;
	    push @gff_files, $ltrc_gff;
	    push @classified_fastas, $ltrc_fas;
	    #unlink $trims_fas, $trims_gff;
	}
    }
    elsif (!$has_trims && $has_ltrs) {
	$log->info("Output files - $trims_gff");
        $log->info("Output files - $trims_fas");
        push @fas_files, $trims_fas;
        push @gff_files, $trims_gff;

	## classifyltrs on LTRs only
        ($ltrc_fas, $ltrc_gff, $ltrc_dir) = $tephra_obj->classify_ltrs($log, $ltr_gff);

	if (-e $ltrc_gff && -e $ltrc_fas) {
            $log->info("Output files - $ltrc_gff");
            $log->info("Output files - $ltrc_fas");
            push @fas_files, $ltrc_fas;
            push @gff_files, $ltrc_gff;
            push @classified_fastas, $ltrc_fas;
            #unlink $trims_fas, $trims_gff;
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

    ## maskref on LTRs
    my $genome_mask2;
    if (defined $ltrc_fas && -e $ltrc_fas && -s $ltrc_fas) {
	$genome_mask2 = $tephra_obj->mask_reference({
            log            => $log,
            config         => $config,
            reference      => $genome_mask1,
            database       => $ltrc_fas,
            dbtype         => 'TRIMs',
            ref_count      => $refct });

        #push @mask_files, $genome_mask2;
	#$refct++;

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
	$log->info("Output files - $hel_gff");
	$log->info("Output files - $hel_fas");
	push @fas_files, $hel_fas;
	#push @gff_files, $hel_gff;
    }

    ## maskref on Helitrons
    my $genome_mask3;
    if (-e $hel_fas && -s $hel_fas) {
	$genome_mask3 = $tephra_obj->mask_reference({
            log            => $log,
            config         => $config,
            reference      => $hel_ref,
            database       => $hel_fas,
            dbtype         => 'Helitrons',
            ref_count      => $refct });

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

    my $tir_gff = $tephra_obj->find_tirs($log, $tir_ref);

    if (-e $tir_gff) {
	$log->info("Output files - $tir_gff");
    }

    ## classifytirs
    my ($tirc_gff, $tirc_fas, $tirc_dir);
    if (-e $tir_gff && -s $tir_gff) {

	($tirc_fas, $tirc_gff) = $tephra_obj->classify_tirs($log, $tir_ref, $tir_gff);

	if (-e $tirc_gff && -e $tirc_fas) {
	    $log->info("Output files - $tirc_gff");
	    $log->info("Output files - $tirc_fas");
	    push @gff_files, $tirc_gff;
	    push @fas_files, $tirc_fas;
	    push @classified_fastas, $tirc_fas;
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
            log            => $log,
            config         => $config,
            reference      => $genome_mask3,
            database       => $tirc_fas,
            dbtype         => 'TIRs',
            ref_count      => $refct });

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
	#push @gff_files, $nonltr_gff;
    }

    ## combine results
    my  ($customRepDB, $customRepGFF) = $tephra_obj->make_combined_repeatdb($log, \@fas_files);
    if (-e $customRepDB) {
        $log->info("Output files - $customRepDB");
    }

    ## maskref on customRepDB
    my $final_mask = $tephra_obj->mask_reference({
            log            => $log,
            config         => $config,
            reference      => $global_opts->{genome},
            database       => $customRepDB,
            dbtype         => 'full transposon database',
            ref_count      => $refct });

    push @mask_files, $final_mask.'.log';

    ## findfragments
    my $fragments_gff = $tephra_obj->find_fragments($log, $customRepDB, $final_mask);

    if (-e $fragments_gff) {
	$log->info("Output files - $fragments_gff");
    }

    ## combine GFF3
    $tephra_obj->combine_gff_files($log, $customRepGFF, \@gff_files, $hel_gff, $nonltr_gff, $fragments_gff);

    ## calculate family similarity
    $tephra_obj->calculate_global_te_similarity($log, $customRepDB);

    ## combine age files
    my $age_ct = $tephra_obj->combine_age_files($log, \@age_files, \@classified_fastas);

    ## clean up
    my $exe_conf = Tephra::Config::Exe->new->get_config_paths;
    my $gt = $exe_conf->{gt};
    my $vmatchbin = $exe_conf->{vmatchbin};
    my $clean_vmidx = File::Spec->catfile($vmatchbin, 'cleanpp.sh');

    $tephra_obj->capture_cmd($clean_vmidx);
    $tephra_obj->capture_cmd($gt, 'clean');
    my @fais = glob "*.fai";
    unlink @fais;
    unlink @mask_files;

    # Log summary of results
    $tephra_obj->log_interval( $tzero, $log );
}

#
# methods
#
sub _get_all_opts {
    my ($config) = @_;

    my ($logfile, $genome, $repeatdb, $hmmdb, $trnadb, $outfile, $clean, $debug, $threads, $subs_rate);
    my ($name, $path, $suffix); # genome file specs
    #dd $config->{all};

    if (defined $config->{all}{genome} && -e $config->{all}{genome}) {
        $genome = $config->{all}{genome};
    }
    else {
	#dd $config->{all};
        say STDERR "\nERROR: genome file was not defined in configuration or does not exist. Check input. Exiting.\n";
        exit(1);
    }

    if (defined $config->{all}{repeatdb} && -e $config->{all}{repeatdb}) {
        $repeatdb = $config->{all}{repeatdb};
    }   
    else {
        say STDERR "\nERROR: repeatdb file was not defined in configuration or does not exist. Check input. Exiting.\n";
        exit(1);
    }

    if (defined $config->{all}{outfile}) {
	$outfile = $config->{all}{outfile};
    }
    else {
	($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
	$outfile = File::Spec->catfile( abs_path($path), $name.'_tephra_transposons.gff3' );
	$logfile = $config->{all}{logfile} // File::Spec->catfile( abs_path($path), $name.'_tephra_full.log' );
    }

    my $execonfig = Tephra::Config::Exe->new->get_config_paths;
    my ($tephra_hmmdb, $tephra_trnadb) = @{$execonfig}{qw(hmmdb trnadb)};
    
    $hmmdb  = defined $config->{all}{hmmdb}  && $config->{all}{hmmdb} =~ /tephradb/i ? 
	$tephra_hmmdb : $config->{all}{hmmdb};
    $trnadb = defined $config->{all}{trnadb} && $config->{all}{trnadb} =~ /tephradb/i ? 
	$tephra_trnadb : $config->{all}{trnadb};

    $clean = defined $config->{all}{clean} && $config->{all}{clean} =~ /yes/i ? 1 : 0;
    $debug = defined $config->{all}{debug} && $config->{all}{debug} =~ /yes/i ? 1 : 0;
    $threads = $config->{all}{threads} // 1;
    $subs_rate = $config->{all}{subs_rate} // 1e-8;
    $logfile = $config->{all}{logfile} // File::Spec->catfile( abs_path($path), $name.'_tephra_full.log' );

    return { logfile   => $logfile,
	     genome    => $genome, 
	     repeatdb  => $repeatdb, 
	     hmmdb     => $hmmdb, 
	     trnadb    => $trnadb, 
	     outfile   => $outfile, 
	     clean     => $clean, 
	     debug     => $debug, 
	     threads   => $threads, 
	     subs_rate => $subs_rate };
}

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
