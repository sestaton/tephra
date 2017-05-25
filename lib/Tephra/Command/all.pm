package Tephra::Command::all;
# ABSTRACT: Run all subcommands and generate annotations for all transposon types.

use 5.014;
use strict;
use warnings;
use File::Spec;
use File::Basename;
use File::Temp;
use Tephra -command;
use Tephra::Config::Reader;
use Tephra::Config::Exe;
use Log::Log4perl;
use Log::Any::Adapter;
use Time::HiRes         qw(gettimeofday);
use POSIX               qw(strftime);
use IPC::System::Simple qw(system capture);
use Cwd                 qw(abs_path);
use Lingua::EN::Inflect;
use DateTime;
use Try::Tiny;
use Cwd;
#use Data::Dump::Color;

our $VERSION = '0.07.3';

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
        $self->help;
        exit(0);
    }
    elsif (! $opt->{config}) {
	say STDERR "\nERROR: Required arguments not given.";
	$self->help and exit(0);
    }
    elsif (! -e $opt->{config}) {
	say STDERR "\nERROR: The configuration file does not exist. Check arguments.";
        $self->help and exit(0);
    }
} 

sub execute {
    my ($self, $opt, $args) = @_;

    my $gff = _run_all_commands($opt);
}

sub _run_all_commands {
    my ($opt) = @_;

    my (@mask_files, @fas_files, @gff_files);
    my $config_obj  = Tephra::Config::Reader->new( config => $opt->{config} );
    my $config      = $config_obj->get_configuration;
    my $global_opts = $config_obj->get_all_opts($config);
    #dd $config and exit;

    ## set global options
    my ($tzero, $log) = _init_tephra($global_opts);

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);   

    ## findltrs
    my $ltr_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs.gff3' );

    my $t0 = gettimeofday();
    my $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findltrs' started at:   $st.");

    my $findltrs_opts = ['-g', $global_opts->{genome}, '-o', $ltr_gff, 
			 '-c', $opt->{config}, '--logfile', $global_opts->{logfile}];
    push @$findltrs_opts, '--debug'
	if $global_opts->{debug};
    _run_tephra_cmd('findltrs', $findltrs_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findltrs' completed at: $ft. Final output file is:");
    $log->info("Output files - $ltr_gff.");

    ## classifyltrs
    my ($ltrc_gff, $ltrc_fas, $ltrc_dir);
    if (-e $ltr_gff && -s $ltr_gff) {
	$ltrc_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs_classified.gff3' );
	$ltrc_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs_classified.fasta' );
	$ltrc_dir = File::Spec->catdir(  abs_path($path), $name.'_tephra_ltrs_classified_results' );

	my $t2 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra classifyltrs' started at:   $st.");
	
	my $classifyltrs_opts = ['-g', $global_opts->{genome}, '-d', $global_opts->{repeatdb}, 
				 '-i', $ltr_gff, '-o', $ltrc_gff,
				 '-r', $ltrc_dir, '-t', $global_opts->{threads}, '--logfile', $global_opts->{logfile}];
	push @$classifyltrs_opts, '--debug'
	    if $global_opts->{debug};
	_run_tephra_cmd('classifyltrs', $classifyltrs_opts, $global_opts->{debug}); 
	
	my $t3 = gettimeofday();
	$total_elapsed = $t3 - $t2;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra classifyltrs' completed at: $ft. Final output files:");
	$log->info("Output files - $ltrc_gff");
	$log->info("Output files - $ltrc_fas.");
	push @fas_files, $ltrc_fas
	    if -e $ltrc_fas;
        push @gff_files, $ltrc_gff
	    if -e $ltrc_gff;
    }

    ## maskref on LTRs
    my $genome_mask1;
    if (-e $ltrc_fas && -s $ltrc_fas) {
	$genome_mask1 = File::Spec->catfile( abs_path($path), $name.'_masked.fasta' );

	my $t4 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' for LTRs started at:   $st.");

	my $mask1_opts = ['-g', $global_opts->{genome}, '-d', $ltrc_fas, '-o', $genome_mask1,
			  '-s', $config->{maskref}{splitsize}, '-v', $config->{maskref}{overlap}, 
			  '-t', $global_opts->{threads}];
	_capture_tephra_cmd('maskref', $mask1_opts, $global_opts->{debug});

	my $t5 = gettimeofday();
	$total_elapsed = $t5 - $t4;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' completed at: $ft. Final output file:");
	$log->info("Output files - $genome_mask1.");
	push @mask_files, $genome_mask1
	    if -e $genome_mask1;
    }

    ## sololtr
    my ($sololtr_gff, $sololtr_rep, $sololtr_fas);
    if (-e $genome_mask1 && -e $ltrc_dir) {
	$sololtr_gff = File::Spec->catfile( abs_path($path), $name.'_sololtrs.gff3' );
	$sololtr_rep = File::Spec->catfile( abs_path($path), $name.'_sololtrs_rep.tsv' );
	$sololtr_fas = File::Spec->catfile( abs_path($path), $name.'_sololtrs_seqs.fasta' );

	my $t6 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra sololtr' started at:   $st.");

	my $solo_opts = ['-i', $ltrc_dir, '-g', $genome_mask1, '-o', $sololtr_gff,
			 '-r', $sololtr_rep, '-s', $sololtr_fas,
			 '-p', $config->{sololtr}{percentid}, '-f', $config->{sololtr}{percentcov},
			 '-l', $config->{sololtr}{matchlen},'-n', $config->{sololtr}{numfamilies},
			 '-t', $global_opts->{threads}];
	push @$solo_opts, '--allfamilies'
	    if $config->{sololtr}{allfamilies} =~ /yes/i;
	_capture_tephra_cmd('sololtr', $solo_opts, $global_opts->{debug});
	
	my $t7 = gettimeofday();
	$total_elapsed = $t7 - $t6;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra sololtr' completed at: $ft.");
	
	if (-e $sololtr_gff && -e $sololtr_rep && -e $sololtr_fas) {
	    $log->info("Final output files:");
	    $log->info("Output files - $sololtr_gff");
	    $log->info("Output files - $sololtr_rep");
	    $log->info("Output files - $sololtr_fas.");
	    push @gff_files, $sololtr_gff;
	}
    }

    ## ltrage
    if (-e $ltrc_gff && -s $ltrc_gff) {
	my $ltrage_out  = File::Spec->catfile( abs_path($path), $name.'_ltrages.tsv' );

	my $t8 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra ltrage' started at:   $st.");

	my $ltrage_opts = ['-g', $global_opts->{genome}, '-t', $global_opts->{threads},
			   '-o', $ltrage_out, '-f', $ltrc_gff, '--clean'];
	push @$ltrage_opts, '--all'
	    if $config->{ltrage}{all} =~ /yes/i;
	if ($config->{ltrage}{all} =~ /no/i) {
	    push @$ltrage_opts, '-i';
	    push @$ltrage_opts, $ltrc_dir;
	}
	_capture_tephra_cmd('ltrage', $ltrage_opts, $global_opts->{debug});
	
	my $t9 = gettimeofday();
	$total_elapsed = $t9 - $t8;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra ltrage' completed at: $ft. Final output file:");
	$log->info("Output files - $ltrage_out.");
    }

    ## illrecomb
    my ($illrec_fas, $illrec_rep, $illrec_stats);
    if (-e $ltrc_fas && -s $ltrc_fas) {
	my $illrec_fas   = File::Spec->catfile( abs_path($path), $name.'_illrecomb.fasta' ); 
	my $illrec_rep   = File::Spec->catfile( abs_path($path), $name.'_illrecomb_rep.tsv' );
	my $illrec_stats = File::Spec->catfile( abs_path($path), $name.'_illrecomb_stats.tsv' );

	my $t10 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra illrecomb' started at:   $st.");

	my $illrec_opts = ['-i', $ltrc_fas, '-o', $illrec_fas, '-r', $illrec_rep,
			   '-s', $illrec_stats, '-t', $global_opts->{threads}];
	_capture_tephra_cmd('illrecomb', $illrec_opts, $global_opts->{debug});
	
	my $t11 = gettimeofday();
	$total_elapsed = $t11 - $t10;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra illrecomb' completed at: $ft. Final output files:");
	$log->info("Output files - $illrec_fas");
	$log->info("Output files - $illrec_rep");
	$log->info("Output files - $illrec_stats.");
    }

    ## TRIMs
    my $trim_ref = (-e $genome_mask1 && -s $genome_mask1) ? $genome_mask1 : $global_opts->{genome};

    my $trims_gff = File::Spec->catfile( abs_path($path), $name.'_trims.gff3' );
    my $trims_fas = File::Spec->catfile( abs_path($path), $name.'_trims.fasta' );
    
    my $t12 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtrims' started at:   $st.");
    
    my $findtrims_opts = ['-g', $trim_ref, '-o', $trims_gff, '--logfile', $global_opts->{logfile}, '--clean'];
    _run_tephra_cmd('findtrims', $findtrims_opts, $global_opts->{debug});
    
    my $t13 = gettimeofday();
    $total_elapsed = $t13 - $t12;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtrims' completed at: $ft. Final output files:");
    $log->info("Output files - $trims_gff");
    $log->info("Output files - $trims_fas.");
    push @fas_files, $trims_fas
	if -e $trims_fas;
    push @gff_files, $trims_gff
	if -e $trims_gff;

    ## maskref for TRIMs
    my $genome_mask2;
    if (-e $trims_fas && -s $trims_fas) {
	$genome_mask2 = File::Spec->catfile( abs_path($path), $name.'_masked2.fasta' );

	my $t14 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' on TRIMs started at:   $st.");

	my $mask2_opts = ['-g', $trim_ref, '-d', $trims_fas, '-o', $genome_mask2,
			  '-s', $config->{maskref}{splitsize}, '-v', $config->{maskref}{overlap},
			  '-t', $global_opts->{threads}];
	_capture_tephra_cmd('maskref', $mask2_opts, $global_opts->{debug});
	
	my $t15 = gettimeofday();
	$total_elapsed = $t15 - $t14;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' on TRIMs completed at: $ft. Final output file:");
	$log->info("Output files - $genome_mask2.");
	push @mask_files, $genome_mask2
	    if -e $genome_mask2;
    }

    ## findhelitrons
    my $hel_ref = (defined $genome_mask2 && -s $genome_mask2) ? $genome_mask2 
	        : (defined $genome_mask1 && -s $genome_mask1) ? $genome_mask1 
	        : $global_opts->{genome};

    my $hel_gff = File::Spec->catfile( abs_path($path), $name.'_helitrons.gff3' );
    my $hel_fas = File::Spec->catfile( abs_path($path), $name.'_helitrons.fasta' );
    
    my $t16 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findhelitrons' started at:   $st.");

    my $findhels_opts = ['-g', $hel_ref, '-o', $hel_gff];
    push @$findhels_opts, '--debug'
	if $global_opts->{debug};
    _capture_tephra_cmd('findhelitrons', $findhels_opts, $global_opts->{debug});

    my $t17 = gettimeofday();
    $total_elapsed = $t17 - $t16;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findhelitrons' completed at: $ft. Final output files:");
    $log->info("Output files - $hel_gff");
    $log->info("Output files - $hel_fas.");
    push @fas_files, $hel_fas
	if -e $hel_fas;

    ## maskref on Helitrons
    my $genome_mask3;
    if (-e $hel_fas && -s $hel_fas) {
	$genome_mask3 = File::Spec->catfile( abs_path($path), $name.'_masked3.fasta' );

	my $t18 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' on Helitrons started at:   $st.");

	my $mask3_opts = ['-g', $hel_ref, '-d', $hel_fas, '-o', $genome_mask3,
			  '-s', $config->{maskref}{splitsize}, '-v', $config->{maskref}{overlap},
			  '-p', $config->{maskref}{percentid}, '-l', $config->{maskref}{hitlength},
			  '-t', $global_opts->{threads}];
	_capture_tephra_cmd('maskref', $mask3_opts, $global_opts->{debug});
	
	my $t19 = gettimeofday();
	$total_elapsed = $t19 - $t18;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' on Helitrons completed at: $ft. Final output file:");
	$log->info("Output files - $genome_mask3.");
	push @mask_files, $genome_mask3
	    if -e $genome_mask3;
    }

    ## findtirs
    my $tir_ref = (defined $genome_mask3 && -s $genome_mask3) ? $genome_mask3
	        : (defined $genome_mask2 && -s $genome_mask2) ? $genome_mask2 
	        : (defined $genome_mask1 && -s $genome_mask1) ? $genome_mask1 
	        : $global_opts->{genome};

    my $t20 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtirs' started at:   $st.");
    my $tir_gff = File::Spec->catfile( abs_path($path), $name.'_tirs.gff3' );

    my $findtirs_opts = ['-g', $tir_ref, '-o', $tir_gff];
    push @$findtirs_opts, '--debug'
        if $global_opts->{debug};
    _capture_tephra_cmd('findtirs', $findtirs_opts, $global_opts->{debug});

    my $t21 = gettimeofday();
    $total_elapsed = $t21 - $t20;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtirs' completed at: $ft. Final output file:");
    $log->info("Output files - $tir_gff.");

    ## classifytirs
    my ($tirc_gff, $tirc_fas);
    if (-e $tir_gff && -s $tir_gff) {
	$tirc_gff = File::Spec->catfile( abs_path($path), $name.'_tirs_classified.gff3' );
	$tirc_fas = File::Spec->catfile( abs_path($path), $name.'_tirs_classified.fasta' );

	my $t22 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra classifytirs' started at:   $st.");

	my $classifytirs_opts = ['-g', $tir_ref, '-i', $tir_gff, '-o', $tirc_gff, '--logfile', $global_opts->{logfile}];
	_run_tephra_cmd('classifytirs', $classifytirs_opts, $global_opts->{debug});

	my $t23 = gettimeofday();
	$total_elapsed = $t23 - $t22;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra classifytirs' completed at: $ft. Final output files:");
	$log->info("Output files - $tirc_gff");
	$log->info("Output files - $tirc_fas.");
	push @gff_files, $tirc_gff
	    if -e $tirc_gff;
        push @fas_files, $tirc_fas
	    if -e $tirc_fas;
    }

    ## maskref on TIRs
    my $genome_mask4;
    if (-e $tirc_fas && -s $tirc_fas) {
	$genome_mask4 = File::Spec->catfile( abs_path($path), $name.'_masked4.fasta' );

	my $t24 = gettimeofday();
	$st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' on TIRs started at:   $st.");
	
	my $mask4_opts = ['-g', $genome_mask3, '-d', $tirc_fas, '-o', $genome_mask4,
			  '-s', $config->{maskref}{splitsize}, '-v', $config->{maskref}{overlap},
			  '-t', $global_opts->{threads}];
	_capture_tephra_cmd('maskref', $mask4_opts, $global_opts->{debug});
	
	my $t25 = gettimeofday();
	$total_elapsed = $t25 - $t24;
	$final_time = sprintf("%.2f",$total_elapsed/60);
	$ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
	$log->info("Command - 'tephra maskref' on TIRs completed at: $ft. Final output file:");
	$log->info("Output files - $genome_mask4.");
	push @mask_files, $genome_mask4
	    if -e $genome_mask4;
    }

    ## findnonltrs
    my $nonltr_ref = (defined $genome_mask4 && -s $genome_mask4) ? $genome_mask4
	           : (defined $genome_mask3 && -s $genome_mask3) ? $genome_mask3
	           : (defined $genome_mask2 && -s $genome_mask2) ? $genome_mask2 
	           : (defined $genome_mask1 && -s $genome_mask1) ? $genome_mask1 
	           : $global_opts->{genome};

    my $t26 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findnonltrs' started at:   $st.");
    my $nonltr_gff = File::Spec->catfile( abs_path($path), $name.'_nonLTRs.gff3' );
    my $nonltr_fas = File::Spec->catfile( abs_path($path), $name.'_nonLTRs.fasta' );

    my $findnonltrs_opts = ['-g', $nonltr_ref, '-o', $nonltr_gff];
    _capture_tephra_cmd('findnonltrs', $findnonltrs_opts, $global_opts->{debug});
    
    my $t27 = gettimeofday();
    $total_elapsed = $t27 - $t26;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findnonltrs' completed at: $ft. Final output files:");
    $log->info("Output files - $nonltr_gff");
    $log->info("Output files - $nonltr_fas.");
    push @fas_files, $nonltr_fas
	if -e $nonltr_fas;
    push @gff_files, $nonltr_gff
	if -e $nonltr_gff;

    ## combine results
    my $t28 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - Generating combined FASTA file at:             $st.");
    my $customRepDB = $global_opts->{outfile} =~ s/\.gff.*/.fasta/r;
    my $tmpiname = $name.'_tephra_transposons_tmp_XXXX';
    my $customRepGFFfh = File::Temp->new( TEMPLATE => $tmpiname,
					  DIR      => abs_path($path),
					  SUFFIX   => '.gff3',
					  UNLINK   => 0);
    my $customRepGFF = $customRepGFFfh->filename;
    open my $out, '>', $customRepDB or die "\nERROR: Could not open file: $customRepDB\n";

    for my $file (@fas_files) {
	my $lines = do { 
	    local $/ = undef; 
	    open my $fh_in, '<', $file or die "\nERROR: Could not open file: $file\n";
	    <$fh_in>;
	};
	print $out $lines;
    }
    close $out;

    my $t29 = gettimeofday();
    $total_elapsed = $t29 - $t28;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - Finished generating combined FASTA file at:     $ft. Final output files:");
    $log->info("Output files - $customRepDB.");
    
    ## maskref on customRepDB
    my $final_mask = File::Spec->catfile( abs_path($path), $name.'_tephra_genome_masked.fasta' );
    my $t30 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra maskref' on full transposon database started at:   $st.");
 
    my $finalmask_opts = ['-g', $global_opts->{genome}, '-d', $customRepDB, '-o', $final_mask,
			  '-s', $config->{maskref}{splitsize}, '-v', $config->{maskref}{overlap},
			  '-t', $global_opts->{threads}];
    _run_tephra_cmd('maskref', $finalmask_opts, $global_opts->{debug});

    my $t31 = gettimeofday();
    $total_elapsed = $t31 - $t30;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - 'tephra maskref' on full transposon database finished at:  $ft. Final output file:");
    $log->info("Output files - $final_mask.");

    ## findfragments
    my $t32 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findfragments' started at:   $st.");
    my $fragments_gff = File::Spec->catfile( abs_path($path), $name.'_transposon_fragments.gff3' );

    my $findfragments_opts = ['-g', $final_mask, '-d', $customRepDB, '-o', $fragments_gff];
    _capture_tephra_cmd('findfragments', $findfragments_opts, $global_opts->{debug});
    
    my $t33 = gettimeofday();
    $total_elapsed = $t33 - $t32;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findfragments' completed at: $ft. Final output files:");
    $log->info("Output files - $fragments_gff");
    push @gff_files, $fragments_gff
	if -e $fragments_gff;

    ## combine GFF3
    my $t34 = gettimeofday();
    $st = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - Generating combined GFF3 file at:              $st.");

    @gff_files = grep { -s $_ } @gff_files; # remove empty files
    my $exe_conf = Tephra::Config::Exe->new->get_config_paths;
    my $gt = $exe_conf->{gt};
    my $gff_cmd = "$gt gff3 -sort @gff_files";
    $gff_cmd .= " | perl -ne 'print unless /^#\\w+\\d+?\$/' > $customRepGFF";
    #say STDERR "debug: $gff_cmd";

    # this is so gt does not drop the Helitron IDs
    my @gtsort_out = capture([0..5], $gff_cmd);
    $gff_cmd = "$gt gff3 -sort -retainids $customRepGFF $hel_gff > $global_opts->{outfile}";
    @gtsort_out = capture([0..5], $gff_cmd);
    unlink $customRepGFF;

    my $t35 = gettimeofday();
    $total_elapsed = $t35 - $t34;
    $final_time = sprintf("%.2f",$total_elapsed/60);
    $ft = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - Finished generating combined GFF3 file at:      $ft. Final output files:");
    $log->info("Output files - $global_opts->{outfile}");

    ## clean up
    my $clean_vmidx;
    my $has_vmatch = 0;
    my @path = split /:|;/, $ENV{PATH};
    for my $p (@path) {
	my $exe = File::Spec->catfile($p, 'cleanpp.sh');

	if (-e $exe && -x $exe) {
	    $clean_vmidx = $exe;
	    $has_vmatch = 1;
	}
    }

    capture([0..5], $clean_vmidx) if $has_vmatch;
    capture([0..5], $gt, 'clean');
    my @fais = glob "*.fai";
    unlink @fais;
    unlink @mask_files;

    # Log summary of results
    _log_interval( $tzero, $log );
}

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

sub _init_tephra {
    my ($config) = @_;

    #my $log_file = File::Spec->catfile($config->{output_directory}, $config->{run_log_file});
    my $conf = qq{
    log4perl.category.Tephra      = INFO, Logfile, Screen
    log4perl.appender.Logfile          = Log::Log4perl::Appender::File
    log4perl.appender.Logfile.filename = $config->{logfile}
    log4perl.appender.Logfile.layout   = Log::Log4perl::Layout::PatternLayout
    log4perl.appender.Logfile.layout.ConversionPattern = %m%n
    log4perl.appender.Screen         = Log::Log4perl::Appender::Screen
    log4perl.appender.Screen.stderr  = 1
    log4perl.appender.Screen.layout  = Log::Log4perl::Layout::SimpleLayout
    };

    Log::Log4perl::init( \$conf );
    Log::Any::Adapter->set('Log4perl');
    my $log = Log::Any->get_logger( category => "Tephra" );
    
    my $t0 = [Time::HiRes::gettimeofday()];
    my $ts = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("======== Tephra version: $VERSION (started at: $ts) ========");
    $log->info("Configuration - Log file for monitoring progress and errors: $config->{logfile}");
    $log->info("Configuration - Genome file:                                 $config->{genome}");
    $log->info("Configuration - Repeat database:                             $config->{repeatdb}");
    $log->info("Configuration - Number of threads:                           $config->{threads}");

    return ($t0, $log);
}

sub _log_interval {
    my ($t0, $log) = @_;
    
    #load_classes('DateTime', 'Time::HiRes', 'Lingua::EN::Inflect', 'POSIX');

    my $t1    = [Time::HiRes::gettimeofday()];
    my $t0_t1 = Time::HiRes::tv_interval($t0, $t1);
    my $dt    = DateTime->from_epoch( epoch => 0 );

    $dt = $dt->add( seconds => $t0_t1 );
    $dt = $dt - DateTime->from_epoch( epoch => 0 );
    
    my @time;
    push @time, $dt->days . Lingua::EN::Inflect::PL_N( ' day', $dt->days ) if $dt->days;
    push @time, $dt->hours . Lingua::EN::Inflect::PL_N( ' hour', $dt->hours ) if $dt->hours;
    push @time, $dt->minutes . Lingua::EN::Inflect::PL_N( ' minute', $dt->minutes ) if $dt->minutes;
    push @time, $dt->seconds . Lingua::EN::Inflect::PL_N( ' second', $dt->seconds ) if $dt->seconds;
    my $timestr = join ', ', @time;
    
    my $fs = POSIX::strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("======== Tephra completed at: $fs. Elapsed time: $timestr. ========");
}

sub _run_tephra_cmd {
    my ($subcmd, $opts, $debug) = @_;

    say STDERR join q{ }, 'DEBUG: ', 'tephra', $subcmd, @$opts 
	if $debug;

    try {
	my @makedbout = system([0..5], 'tephra', $subcmd, @$opts);
    }
    catch {
	say STDERR "Unable to run 'tephra $subcmd'. Here is the exception: $_.";
	exit(1);
    };
}

sub _capture_tephra_cmd {
    my ($subcmd, $opts, $debug) = @_;

    say STDERR join q{ }, 'DEBUG: ', 'tephra', $subcmd, @$opts
	if $debug;

    try {
	my @makedbout = capture([0..5], 'tephra', $subcmd, @$opts);
    }
    catch {
	say STDERR "Unable to run 'tephra $subcmd'. Here is the exception: $_.";
	exit(1);
    };
}
   
sub help {
    print STDERR<<END

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
 along with a FASTA database and masked genome. An table of results showing global repeat composition 
 will be produced for the input genome.

=head1 AUTHOR 

S. Evan Staton, C<< <statonse at gmail.com> >>

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
