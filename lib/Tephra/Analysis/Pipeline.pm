package Tephra::Analysis::Pipeline;

use 5.014;
use Moose;
use Cwd         qw(abs_path);
use Time::HiRes qw(gettimeofday);
use POSIX       qw(strftime);
use File::Temp  qw(tempfile);
use File::Basename;
use File::Spec;
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use Tephra::Config::Exe;
use Tephra::Annotation::Util;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::Logger',
     'Tephra::Role::Run::Any',
     'Tephra::Role::Run::TephraCmd';

=head1 NAME

Tephra::Analysis::Pipeline - Methods for running the 'tephra all' command

=head1 VERSION

Version 0.09.9

=cut

our $VERSION = '0.09.9';
$VERSION = eval $VERSION;

has global_options => ( is => 'ro', isa => 'HashRef', required => 1 );
has cli_options => ( is => 'ro', required => 1 );

sub find_ltrs {
    my $self = shift;
    my ($log) = @_;
    my $global_opts = $self->global_options;
    my $opt = $self->cli_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/); 
    my $ltr_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs.gff3' );
    my $ltr_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs.fasta' );

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findltrs' started at: $st.");

    my $findltrs_opts = ['-g', $global_opts->{genome}, '-o', $ltr_gff, 
                         '-c', $opt->{config}, '--logfile', $global_opts->{logfile}];
    push @$findltrs_opts, '--debug'
        if $global_opts->{debug};
    $self->run_tephra_cmd('findltrs', $findltrs_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findltrs' completed at: $ft.");

    return ($ltr_fas, $ltr_gff);
}

sub find_trims {
    my $self = shift;
    my ($log, $trim_ref) = @_;
    my $global_opts = $self->global_options;
    my $opt = $self->cli_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/); 
    my $trims_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_trims.gff3' );
    my $trims_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_trims.fasta' );
    
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtrims' started at: $st.");
    
    my $findtrims_opts = ['-g', $trim_ref, '-o', $trims_gff, '--logfile', $global_opts->{logfile}, '--clean'];
    $self->run_tephra_cmd('findtrims', $findtrims_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtrims' completed at: $ft.");

    return ($trims_fas, $trims_gff);
}

sub mask_reference {
    my $self = shift;
    my ($data_ref) = @_;
    my $global_opts = $self->global_options;

    my ($log, $log_results, $config, $ref, $db, $type, $refct) = 
        @{$data_ref}{qw(log log_results config reference database dbtype ref_count)};

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/); 

    my $file;
    if ($type =~ /full/) {
        $file = $name.'_tephra_genome_masked.fasta';
    }
    else {
        if ($refct > 1) {
            $file = $name."_tephra_masked$refct.fasta";
        }
        else {
            $file = $name.'_tephra_masked.fasta';
        }
    }

    my $masked_ref = File::Spec->catfile( abs_path($path), $file );

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra maskref' for $type started at: $st.");
    
    my $mask_opts = ['-g', $ref, '-d', $db, '-o', $masked_ref, '-s', $config->{maskref}{splitsize}, 
		     '-v', $config->{maskref}{overlap}, '-p', $config->{maskref}{percentid}, 
		     '-l', $config->{maskref}{hitlength}, '-t', $global_opts->{threads}];

    if ($log_results) {
	$self->run_tephra_cmd('maskref', $mask_opts, $global_opts->{debug});
    }
    else {
	$self->capture_tephra_cmd('maskref', $mask_opts, $global_opts->{debug});
    }

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra maskref' completed at: $ft. Final output file:");
    
    if (-e $masked_ref) { 
        return $masked_ref;
    }
    else {
        return undef;
    }
}

sub combine_ltrs_trims {
    my $self = shift;
    my ($trims_gff, $ltr_gff) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $ltr_trim_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs_trims.fasta' );
    my $exe_conf = Tephra::Config::Exe->new->get_config_paths;
    my $gt = $exe_conf->{gt};
    my $gff_cmd = "$gt gff3 -sort $trims_gff $ltr_gff > $ltr_trim_gff";

    $self->capture_cmd($gff_cmd);

    return $ltr_trim_gff;
}

sub classify_ltrs {
    my $self = shift;
    my ($log, $ltr_trim_gff) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $ltrc_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs_trims_classified.gff3' );
    my $ltrc_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrs_trims_classified.fasta' );
    my $ltrc_dir = File::Spec->catdir(  abs_path($path), $name.'_tephra_ltrs_trims_classified_results' );

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra classifyltrs' started at: $st.");
        
    my $classifyltrs_opts = ['-g', $global_opts->{genome}, '-d', $global_opts->{repeatdb}, '-i', $ltr_trim_gff, 
                             '-o', $ltrc_gff, '-r', $ltrc_dir, '-t', $global_opts->{threads}, 
                             '--logfile', $global_opts->{logfile}];
    push @$classifyltrs_opts, '--debug'
	if $global_opts->{debug};
    $self->run_tephra_cmd('classifyltrs', $classifyltrs_opts, $global_opts->{debug}); 
        
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra classifyltrs' completed at: $ft.");

    return ($ltrc_fas, $ltrc_gff, $ltrc_dir);
}

sub calculate_ltrage {
    my $self = shift;
    my ($config, $log, $ltrc_gff, $ltrc_dir) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $ltrage_out = File::Spec->catfile( abs_path($path), $name.'_tephra_ltrages.tsv' );

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra ltrage' started at: $st.");

    my $ltrage_opts = ['-g', $global_opts->{genome}, '-t', $global_opts->{threads},
                       '-o', $ltrage_out, '-f', $ltrc_gff, '-r', $global_opts->{subs_rate}, '--clean'];
        push @$ltrage_opts, '--all'
            if $config->{ltrage}{all} =~ /yes/i;

    if ($config->{ltrage}{all} =~ /no/i) {
        push @$ltrage_opts, '-i';
        push @$ltrage_opts, $ltrc_dir;
    }
    $self->capture_tephra_cmd('ltrage', $ltrage_opts, $global_opts->{debug});
        
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra ltrage' completed at: $ft.");

    return $ltrage_out;
}

sub find_sololtrs {
    my $self = shift;
    my ($config, $log, $genome_mask1, $ltrc_dir) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $sololtr_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_sololtrs.gff3' );
    my $sololtr_rep = File::Spec->catfile( abs_path($path), $name.'_tephra_sololtrs_rep.tsv' );
    my $sololtr_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_sololtrs_seqs.fasta' );

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra sololtr' started at: $st.");

    my $solo_opts = ['-i', $ltrc_dir, '-g', $genome_mask1, '-o', $sololtr_gff,
                     '-r', $sololtr_rep, '-s', $sololtr_fas,
                     '-p', $config->{sololtr}{percentid}, '-f', $config->{sololtr}{percentcov},
                     '-l', $config->{sololtr}{matchlen},'-n', $config->{sololtr}{numfamilies},
                     '-t', $global_opts->{threads}];
    push @$solo_opts, '--allfamilies'
        if $config->{sololtr}{allfamilies} =~ /yes/i;
    $self->capture_tephra_cmd('sololtr', $solo_opts, $global_opts->{debug});
        
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra sololtr' completed at: $ft.");
    
    return ($sololtr_gff, $sololtr_rep, $sololtr_fas);
}

sub find_illrecombination {
    my $self = shift;
    my ($log, $ltrc_fas) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $illrec_fas   = File::Spec->catfile( abs_path($path), $name.'_tephra_illrecomb.fasta' ); 
    my $illrec_rep   = File::Spec->catfile( abs_path($path), $name.'_tephra_illrecomb_rep.tsv' );
    my $illrec_stats = File::Spec->catfile( abs_path($path), $name.'_tephra_illrecomb_stats.tsv' );
    
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra illrecomb' started at: $st.");
    
    my $illrec_opts = ['-i', $ltrc_fas, '-o', $illrec_fas, '-r', $illrec_rep,
                       '-s', $illrec_stats, '-t', $global_opts->{threads}];
    $self->capture_tephra_cmd('illrecomb', $illrec_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra illrecomb' completed at: $ft.");

    return ($illrec_fas, $illrec_rep, $illrec_stats);
}

sub find_helitrons {
    my $self = shift;
    my ($log, $hel_ref) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $hel_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_helitrons.gff3' );
    my $hel_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_helitrons.fasta' );
    
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findhelitrons' started at: $st.");

    my $findhels_opts = ['-g', $hel_ref, '-o', $hel_gff, '--logfile', $global_opts->{logfile}];
    push @$findhels_opts, '--debug'
        if $global_opts->{debug};
    $self->capture_tephra_cmd('findhelitrons', $findhels_opts, $global_opts->{debug});

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findhelitrons' completed at: $ft.");

    return ($hel_fas, $hel_gff);
}

sub find_tirs {
    my $self = shift;
    my ($log, $tir_ref) = @_;
    my $global_opts = $self->global_options;

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtirs' started at: $st.");

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $tir_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_tirs.gff3' );
    #my $tir_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_tirs.fasta' );

    my $findtirs_opts = ['-g', $tir_ref, '-o', $tir_gff];
    push @$findtirs_opts, '--debug'
        if $global_opts->{debug};
    $self->capture_tephra_cmd('findtirs', $findtirs_opts, $global_opts->{debug});

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findtirs' completed at: $ft.");

    return $tir_gff;
}

sub classify_tirs {
    my $self = shift;
    my ($log, $tir_ref, $tir_gff) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $tirc_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_tirs_classified.gff3' );
    my $tirc_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_tirs_classified.fasta' );
    my $tirc_dir = File::Spec->catdir(  abs_path($path), $name.'_tephra_tirs_classified_results' );
    
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra classifytirs' started at: $st.");

    my $classifytirs_opts = ['-g', $global_opts->{genome}, '-i', $tir_gff, '-o', $tirc_gff, '-r', $tirc_dir, 
                             '-t', $global_opts->{threads}, '-d', $global_opts->{repeatdb}, 
                             '--logfile', $global_opts->{logfile}];
    $self->run_tephra_cmd('classifytirs', $classifytirs_opts, $global_opts->{debug});

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra classifytirs' completed at: $ft.");

    return ($tirc_fas, $tirc_gff, $tirc_dir);
}

sub calculate_tirage {
    my $self = shift;
    my ($config, $log, $tirc_gff, $tirc_dir) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $tirage_out = File::Spec->catfile( abs_path($path), $name.'_tephra_tirages.tsv' );
    
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra tirage' started at: $st.");
    
    my $tirage_opts = ['-g', $global_opts->{genome}, '-t', $global_opts->{threads},
                       '-o', $tirage_out, '-f', $tirc_gff, '-r', $global_opts->{subs_rate}, '--clean'];
    push @$tirage_opts, '--all'
        if $config->{tirage}{all} =~ /yes/i;
    if ($config->{tirage}{all} =~ /no/i) {
        push @$tirage_opts, '-i';
        push @$tirage_opts, $tirc_dir;
    }
    $self->capture_tephra_cmd('tirage', $tirage_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra tirage' completed at: $ft.");
    
    return $tirage_out;
}

sub find_nonltrs {
    my $self = shift;
    my ($log, $nonltr_ref) = @_;
    my $global_opts = $self->global_options;

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findnonltrs' started at: $st.");

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $nonltr_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_nonLTRs.gff3' );
    my $nonltr_fas = File::Spec->catfile( abs_path($path), $name.'_tephra_nonLTRs.fasta' );

    my $findnonltrs_opts = ['-g', $nonltr_ref, '-r', $global_opts->{genome}, '-o', $nonltr_gff, 
			    '--logfile', $global_opts->{logfile}];
    $self->capture_tephra_cmd('findnonltrs', $findnonltrs_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findnonltrs' completed at: $ft.");

    return ($nonltr_fas, $nonltr_gff);
}

sub make_combined_repeatdb {
    my $self = shift;
    my ($log, $fas_files) = @_;
    my $global_opts = $self->global_options;

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - Generating combined FASTA file at: $st.");
    
    my $customRepDB = $global_opts->{outfile} =~ s/\.gff.*/.fasta/r;
    open my $out, '>', $customRepDB or die "\n[ERROR]: Could not open file: $customRepDB\n";

    for my $file (@$fas_files) {
        my $lines = do { 
            local $/ = undef; 
            open my $fh_in, '<', $file or die "\n[ERROR]: Could not open file: $file\n";
            <$fh_in>;
        };
        print $out $lines;
    }
    close $out;

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - Finished generating combined FASTA file at: $ft. Final output files:");

    return $customRepDB;
}

sub find_fragments {
    my $self = shift;
    my ($log, $customRepDB, $final_mask) = @_;
    my $global_opts = $self->global_options;

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findfragments' started at: $st.");

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $fragments_gff = File::Spec->catfile( abs_path($path), $name.'_tephra_transposon_fragments.gff3' );

    my $findfragments_opts = ['-g', $final_mask, '-d', $customRepDB, '-o', $fragments_gff, '-t', $global_opts->{threads}];
    $self->capture_tephra_cmd('findfragments', $findfragments_opts, $global_opts->{debug});
    
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - 'tephra findfragments' completed at: $ft.");

    return $fragments_gff;
}

sub combine_gff_files {
    my $self = shift;
    my ($log, $gff_files, $hel_gff, $nonltr_gff, $fragments_gff) = @_;
    my $global_opts = $self->global_options;

    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - Generating combined GFF3 file at: $st.");

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $tmpiname = $name.'_tephra_transposons_tmp_XXXX';
    my ($outfh, $customRepGFF) = tempfile( TEMPLATE => $tmpiname, DIR => abs_path($path), SUFFIX => '.gff3', UNLINK => 0 );
    close $outfh;
    @$gff_files = grep { -s $_ } @$gff_files; # remove empty files
    my $exe_conf = Tephra::Config::Exe->new->get_config_paths;
    my $gt = $exe_conf->{gt};
    my $gff_cmd = "$gt gff3 -sort -retainids @$gff_files";
    $gff_cmd .= " | perl -ne 'print unless /^#\\w+\\d+?\$/' > $customRepGFF";
    $self->capture_cmd($gff_cmd);

    # this is so gt does not drop IDs
    my @final_gffs = grep { -s $_ } ($customRepGFF, $hel_gff, $nonltr_gff, $fragments_gff); # remove empty files
    $gff_cmd = "$gt gff3 -sort -retainids @final_gffs > $global_opts->{outfile}";
    say STDERR "debug: $gff_cmd";
    $self->capture_cmd($gff_cmd);
    unlink $customRepGFF;

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - Finished generating combined GFF3 file at: $ft. Final output files:");

    if (-e $global_opts->{outfile}) {
        $log->info("Output files - $global_opts->{outfile}");
    }

    return;
}

sub calculate_global_te_similarity {
    my $self = shift;
    my ($log, $customRepDB) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $te_sum = File::Spec->catfile( abs_path($path), $name.'_tephra_transposons_length-similarity_summary.tsv' );
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Command - Calculating global transposon family similarity at: $st.");

    my $util = Tephra::Annotation::Util->new;
    $util->calculate_family_similarity($customRepDB, $te_sum, $global_opts->{threads});

    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - Finished calculating global transposon family similarity at: $st.");

    return;
}

sub combine_age_files {
    my $self = shift;
    my ($log, $ages, $fastas) = @_;
    my $global_opts = $self->global_options;

    my ($name, $path, $suffix) = fileparse($global_opts->{genome}, qr/\.[^.]*/);
    my $age_sum = File::Spec->catfile( abs_path($path), $name.'_tephra_ltr-tir_age_summary.tsv' );
    my $t0 = gettimeofday();
    my $st = strftime('%d-%m-%Y %H:%M:%S', localtime);

    my $util = Tephra::Annotation::Util->new;

    my ($famname, %families, %ages);
    for my $fasta (@$fastas) {
        open my $in, '<', $fasta or die "\n[ERROR]: Could not open file: $fasta\n";
        
        while (my $line = <$in>) {
            chomp $line;
            if ($line =~ /^>([A-Z]{3}(?:_singleton_)?(?:_?family\d+)?)_(?:LTR|TRIM|terminal)/) {
                $famname = $1;
                $line =~ s/>//;
                push @{$families{$famname}}, $line;
            }
        }
        close $in;
    }

    for my $agefile (@$ages) { 
        open my $tab, '<', $agefile or die "\n[ERROR]: Could not open file: $agefile\n";
        
        while (my $line = <$tab>) {
            chomp $line;
            next if $line =~ /^(?:LTR|TIR)-ID/;
            #LTR-ID Divergence Age Ts:Tv
            my ($id, $div, $age, $tstv) = split /\t/, $line;
            if ($id =~ /^([A-Z]{3}(?:_singleton_)?(?:_?family\d+)?)_(?:LTR|TRIM|terminal)/) {
                my $fam = $1;
                if (exists $families{$fam}) {
                    my $famsize = @{$families{$fam}};
                    push @{$ages{$famsize}{$fam}}, join "||", $id, $div, $age, $tstv;
                }
            }
        }
        close $tab;
    }

    open my $out, '>', $age_sum or die "\n[ERROR]: Could not open file: $age_sum\n";
    say $out join "\t", 'Superfamily', 'Family', 'Family_size', 'ElementID', 'Divergence', 'Age', 'Ts:Tv';

    my $ct = 0;
    for my $famsize (reverse sort { $a <=> $b } keys %ages) {
        for my $fam (nsort keys %{$ages{$famsize}}) { 
            for my $agestr (@{$ages{$famsize}{$fam}}) {
                $ct++;
                my ($id, $div, $age, $tstv) = split /\|\|/, $agestr;
                my $sfam = $util->map_superfamily_name($fam);
                $sfam = $sfam ? $sfam : 'Unknown'; # the method above returns 0 if the superfamily is unknown
                say $out join "\t", $sfam, $fam, $famsize, $id, $div, $age, $tstv;
            }
        }
    }
    close $out;

    my $age_pd = ' ' x length($ct);
    $log->info("Command - Generating combined age files at: $st.");
    my $t1 = gettimeofday();
    my $total_elapsed = $t1 - $t0;
    my $final_time = sprintf("%.2f",$total_elapsed/60);
    my $ft = strftime('%d-%m-%Y %H:%M:%S', localtime);
    $log->info("Results - Finished generating combined age files for $ct transposons at: $st.");

    return $ct;
}

sub make_temp_reference_for_masking {
    my $self = shift;
    my ($ltr_fas, $te_type) = @_;

    # LTR or TRIM
    my $code = $te_type =~ /ltr/i ? 'RLX' : 'RLT';

    my ($name, $path, $suffix) = fileparse($ltr_fas, qr/\.[^.]*/);
    my $tmpiname = $name.'_tephra_'.$te_type.'_tmp_XXXX';

    my ($outfh, $tmp_file) = tempfile( TEMPLATE => $tmpiname, DIR => abs_path($path), SUFFIX => '.fasta', UNLINK => 0 );

    my $kseq = Bio::DB::HTS::Kseq->new($ltr_fas);
    my $iter = $kseq->iterator;

    while (my $seqobj = $iter->next_seq) { 
        my $id = join "_", $code, $seqobj->name;
        say $outfh join "\n", ">".$id, $seqobj->seq;
    }
    close $outfh;

    return $tmp_file;
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

    perldoc Tephra::Analysis::Pipeline


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
