package Tephra::Annotation::Util;

use 5.014;
use Moose;
use List::Util qw(sum);
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use namespace::autoclean;
#use Data::Dump::Color;

with 'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Annotation::Util - Utility methods for working with transposon annotations

=head1 VERSION

Version 0.12.5

=cut

our $VERSION = '0.12.5';
$VERSION = eval $VERSION;

has 'debug' => (
    is         => 'ro',
    isa        => 'Bool',
    predicate  => 'has_debug',
    lazy       => 1,
    default    => 0,
);

sub calculate_family_similarity {
    my $self = shift;
    my ($fasta, $outfile, $threads) = @_;

    my $kseq = Bio::DB::HTS::Kseq->new($fasta);
    my $iter = $kseq->iterator;
    
    my ($famname, %families, %lengths);
    while (my $obj = $iter->next_seq) {
	my $id = $obj->name;
	my $seq = $obj->seq;
	#next if $id =~ /TRIM_retrotransposon/i;
	if ($id =~ /^(\w{3}_(?:singleton_)?family\d+)_/) {
	    $famname = $1;
	    push @{$families{$famname}}, { seq => $seq, id => $id };
	    push @{$lengths{$famname}}, length($seq);
	}
	else {
	    say STDERR "\n[WARNING]: '$id' does not match pattern.\n";
	}
    }

    my @results;
    for my $fam (reverse sort { @{$families{$a}} <=> @{$families{$b}} } keys %families) {
	my $famsize = @{$families{$fam}};
	if (@{$families{$fam}} > 1) {
	    my $fas = $fam.'.fasta';
	    open my $out, '>', $fas or die "\n[ERROR]: Could not open file: $fas\n";
	    for my $seqobj (@{$families{$fam}}) {
		say $out join "\n", ">".$seqobj->{id}, $seqobj->{seq};
	    }
	    close $out;
	    my $result = $self->_do_blast_search($fam, $famsize, $fas, $threads, $lengths{$fam});
	    push @results, $result 
		if defined $result;
	    unlink $fas;
	}
    }

    if (@results) {
	open my $outfh, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	say $outfh join "\t", 'family', 'family_size', 'family_count_analyzed', 'similarity_mean', 'hit_length_mean', 'element_length_mean';
	for my $result (@results) {
	    say $outfh $result;
	}
	close $outfh;
    }
    else {
	say STDERR "\n[WARNING]: No multi-member families discovered, so no family-level similarity will be reported.\n";
	unlink $outfile if -e $outfile;
    }

    return;
}

=head2 map_superfamily_name

 Title   : map_superfamily_name

 Usage   : my $superfamily_name = $self->map_superfamily_name($match);
           
 Function: Get the superfamily common name based on the 3-letter code
                                                                            Return_type
 Returns : The TE superfamily name                                          Scalar
            
                                                                            Arg_type
 Args    : A BLAST hit or any sequence identifier                           Scalar

=cut

## NB: This method is pulled straight from Transposome (github.com/sestaton/Transposome).
sub map_superfamily_name {
    my $self = shift;
    my ($id) = @_;

    my ($sfamily_code) = ($id =~ /(^[A-Z]{3})_?/);
    unless (defined $sfamily_code) {
	say STDERR "\n[WARNING]: Could not get 3-letter code from: $id. Skipping.\n"
	    if $self->debug;
	return 0;
    }

    my %sfcode_table = (
        'RYD' => 'DIRS',
        'RYN' => 'Ngaro',
        'RYX' => 'Unknown_DIRS',
        'RYV' => 'VIPER',
	'RIC' => 'CR1', # CR1 clade
        'RII' => 'I',
        'RIJ' => 'Jockey',
        'RIL' => 'L1',
        'RIM' => 'L2',
        'RIS' => 'R1',
        'RIA' => 'RandI',
        'RIC' => 'Rex', # CR1 clade, http://link.springer.com/article/10.1186/1471-2148-13-152/fulltext.html?view=classic
        'RIT' => 'RTE',
        'RID' => 'Tad1',
        'RIR' => 'R2',
        'RIE' => 'CRE', # http://www.ncbi.nlm.nih.gov/pubmed/15939396
        'RIX' => 'Unknown_LINE',
        'RLB' => 'Bel/Pao',
        'RLG' => 'Gypsy',
        'RLC' => 'Copia',
        'RLE' => 'ERV',
        'RLR' => 'Retrovirus',
	'RLT' => 'TRIM',
        'RLX' => 'Unknown_LTR',
        'PPP' => 'Penelope',
        'RPX' => 'Unknown_PLE',
        'RSS' => '5S',
        'RSL' => '7SL',
        'RST' => 'tRNA',
        'RSX' => 'Unknown_SINE',
        'RXX' => 'Unknown_retrotransposon',
        'DYC' => 'Crypton',
        'DYX' => 'Unknown_Crypton',
        'DTC' => 'CACTA',
        'DTA' => 'hAT',
        'DTE' => 'Merlin',
        'DTM' => 'Mutator',
        'DTP' => 'P',
        'DTH' => 'PIF/Harbinger',
        'DTB' => 'PiggyBac',
        'DTT' => 'Tc1/Mariner',
        'DTR' => 'Transib',
	#'DTI' => 'Unknown_MITE',
        'DTX' => 'Unknown_TIR',
        'DXX' => 'Unknown_DNA_transposon',
        'DHH' => 'Helitron',
        'DHX' => 'Unknown_Helitron',
        'DMM' => 'Maverick',
        'DMX' => 'Unknown_Maverick',
	'RST' => 'SINE2/tRNA' );
    
    if (exists $sfcode_table{$sfamily_code}) {
	return $sfcode_table{$sfamily_code};
    }
    else {
	say STDERR "\n[WARNING]: No 3-letter code could be found for: $sfamily_code\n";
	return 0;
    }
}

=head2 map_superfamily_name_to_code

 Title   : map_superfamily_name_to_code

 Usage   : my $superfamily_name = $self->map_superfamily_name_to_code($match);
           
 Function: Get the 3-letter code based on a superfamily or clade name
                                                                            Return_type
 Returns : The 3-letter superfamily code                                    Scalar
            
                                                                            Arg_type
 Args    : A BLAST hit or any sequence identifier                           Scalar

=cut

sub map_superfamily_name_to_code {
    my $self = shift;
    my ($name) = @_;

    my %sfcode_table = (
        'RYD' => 'DIRS',
        'RYN' => 'Ngaro',
        'RYX' => 'Unknown_DIRS',
        'RYV' => 'VIPER',
        'RIX' => 'Unknown_LINE',
        'RLB' => 'Bel/Pao',
        'RLG' => 'Gypsy',
        'RLC' => 'Copia',
        'RLE' => 'ERV',
        'RLR' => 'Retrovirus',
	'RLT' => 'TRIM',
        'RLX' => 'Unknown_LTR',
        'PPP' => 'Penelope',
        'RPX' => 'Unknown_PLE',
        'RSS' => '5S',
        'RSL' => '7SL',
        'RST' => 'tRNA',
        'RSX' => 'Unknown_SINE',
        'RXX' => 'Unknown_retrotransposon',
        'DYC' => 'Crypton',
        'DYX' => 'Unknown_Crypton',
        'DTC' => 'CACTA',
        'DTA' => 'hAT',
        'DTE' => 'Merlin',
        'DTM' => 'Mutator',
        'DTP' => 'P',
        'DTH' => 'PIF/Harbinger',
        'DTB' => 'PiggyBac',
        'DTT' => 'Tc1/Mariner',
        'DTR' => 'Transib',
	#'DTI' => 'Unknown_MITE',
        'DTX' => 'Unknown_TIR',
        'DXX' => 'Unknown_DNA_transposon',
        'DHH' => 'Helitron',
        'DHX' => 'Unknown_Helitron',
        'DMM' => 'Maverick',
        'DMX' => 'Unknown_Maverick',
	'RST' => 'SINE2/tRNA',
        'RIC' => 'CR1', # CR1 clade
        'RII' => 'I',
        'RIJ' => 'Jockey',
        'RIL' => 'L1',
        'RIM' => 'L2',
        'RIS' => 'R1',
        'RIA' => 'RandI',
        'RIC' => 'Rex', # CR1 clade, http://link.springer.com/article/10.1186/1471-2148-13-152/fulltext.html?view=classic      
        'RIT' => 'RTE',
        'RID' => 'Tad1',
        'RIR' => 'R2',
        'RIE' => 'CRE', # http://www.ncbi.nlm.nih.gov/pubmed/15939396
        'RIX' => 'Unknown_LINE',
    );

    my %name_table = reverse %sfcode_table;

    if (exists $name_table{$name}) {
	return $name_table{$name};
    }
    else {
	say STDERR "\n[WARNING]: No 3-letter code could be found for: $name\n";
	return 0;
    }
}

sub build_repeat_map {
    my $self = shift;

    my %repeat_map = (
	## Class I
	# DIRS
	'RLD' => { class => 'Class I', order => 'DIRS', repeat_name => 'DIRS' },
	'RYN' => { class => 'Class I', order => 'DIRS', repeat_name => 'Ngaro' },
	'RYX' => { class => 'Class I', order => 'DIRS', repeat_name => 'Unknown DIRS' },
	'RYV' => { class => 'Class I', order => 'DIRS', repeat_name => 'VIPER' },
	# LINE 
	'RII' => { class => 'Class I', order => 'LINE', repeat_name => 'I' },
	'RIJ' => { class => 'Class I', order => 'LINE', repeat_name => 'Jockey' },
	'RIL' => { class => 'Class I', order => 'LINE', repeat_name => 'L1' },
	'RIR' => { class => 'Class I', order => 'LINE', repeat_name => 'R2' },
	'RIT' => { class => 'Class I', order => 'LINE', repeat_name => 'RTE' },
	'RIX' => { class => 'Class I', order => 'LINE', repeat_name => 'Unknown LINE' },
	'RIC' => { class => 'Class I', order => 'LINE', repeat_name => 'CR1' },
	'RIS' => { class => 'Class I', order => 'LINE', repeat_name => 'R1' },
	'RIA' => { class => 'Class I', order => 'LINE', repeat_name => 'RandI' },
	'RIE' => { class => 'Class I', order => 'LINE', repeat_name => 'Rex' },
	'RID' => { class => 'Class I', order => 'LINE', repeat_name => 'Tad1' },
	# LTR
	'RLB' => { class => 'Class I', order => 'LTR', repeat_name => 'Bel/Pao' },
	'RLC' => { class => 'Class I', order => 'LTR', repeat_name => 'Copia' },
	'RLE' => { class => 'Class I', order => 'LTR', repeat_name => 'ERV' },
	'RLG' => { class => 'Class I', order => 'LTR', repeat_name => 'Gypsy' },
	'RLR' => { class => 'Class I', order => 'LTR', repeat_name => 'Retrovirus' },
	'RLT' => { class => 'Class I', order => 'LTR', repeat_name => 'TRIM' },
	'RLX' => { class => 'Class I', order => 'LTR', repeat_name => 'Unknown LTR' },
	# PLE
	'RPP' => { class => 'Class I', order => 'Penelope', repeat_name => 'Penelope' },
	'RPX' => { class => 'Class I', order => 'Penelope', repeat_name => 'Unknown PLE' },
	# SINE
	'RSS' => { class => 'Class I', order => 'SINE', repeat_name => '5S' },
	'RSL' => { class => 'Class I', order => 'SINE', repeat_name => '7SL' },
	'RST' => { class => 'Class I', order => 'SINE', repeat_name => 'tRNA' },
	'RSX' => { class => 'Class I', order => 'SINE', repeat_name => 'Unknown SINE' },
	'RXX' => { class => 'Class I', order => 'SINE', repeat_name => 'Unknown retrotransposon' },
	## Class II
	# - Subclass 1
	# Crypton
	'DYC' => { class => 'Class II', order => 'Crypton', repeat_name => 'Crypton' },
	'DYX' => { class => 'Class II', order => 'Crypton', repeat_name => 'Unknown Crypton' },
	# TIR
	'DTC' => { class => 'Class II', order => 'TIR', repeat_name => 'CACTA' },
	'DTA' => { class => 'Class II', order => 'TIR', repeat_name => 'hAT' },
	'DTE' => { class => 'Class II', order => 'TIR', repeat_name => 'Merlin' },
	'DTM' => { class => 'Class II', order => 'TIR', repeat_name => 'Mutator' },
	'DTP' => { class => 'Class II', order => 'TIR', repeat_name => 'P' },
	'DTH' => { class => 'Class II', order => 'TIR', repeat_name => 'PIF/Harbinger' },
	'DTB' => { class => 'Class II', order => 'TIR', repeat_name => 'PiggyBac' },
	'DTT' => { class => 'Class II', order => 'TIR', repeat_name => 'Tc1/Mariner' },
	'DTR' => { class => 'Class II', order => 'TIR', repeat_name => 'Transib' },
	#'DTI' => { class => 'Class II', order => 'TIR', repeat_name => 'Unknown MITE' },
	'DTX' => { class => 'Class II', order => 'TIR', repeat_name => 'Unknown TIR' },
	'DXX' => { class => 'Class II', order => 'TIR', repeat_name => 'Unknown DNA transposon' },
	# - Subclass 2
	# Helitron
	'DHH' => { class => 'Class II', order => 'Helitron', repeat_name => 'Helitron' },
	'DHX' => { class => 'Class II', order => 'Helitron', repeat_name => 'Unknown Helitron' },
	# Maverick
	'DMM' => { class => 'Class II', order => 'Maverick', repeat_name => 'Maverick' },
	'DMX' => { class => 'Class II', order => 'Maverick', repeat_name => 'Unknown Maverick' },
	);

    return \%repeat_map;
}

sub get_GO_terms {
    my $self = shift;
    my ($term) = @_;

    my %table = ( 
	'RVT_1' => 'GO:0003964'
    );

    my $has_terms = 0;
    if (exists $table{$term}) {
	$has_terms = 1;

	return $table{$term};
    }
    else {
	say STDERR "\n[WARNING]: No GO terms are defined for '$term' so they will not be added to the ".
	    "GFF3 output.\n";

	return 0;
    }

    return;
}

sub get_SO_terms {
    my $self = shift;

    my %table = (
        'LTR_retrotransposon'     => 'SO:0000186',
        'non_LTR_retrotransposon' => 'SO:0000189',
        
        'U_box'                => 'SO:0001788',
        'RR_tract'             => 'SO:0000435',
        'long_terminal_repeat' => 'SO:0000286',
        'inverted_repeat'      => 'SO:0000294',
        'primer_binding_site'  => 'SO:0005850',
        'protein_match'        => 'SO:0000349',
        
        'terminal_inverted_repeat_element' => 'SO:0000208',
        'terminal_inverted_repeat'         => 'SO:0000481',
        'helitron'                         => 'SO:0000544',
        'MITE'                             => 'SO:0000338',
        'DNA_transposon'                   => 'SO:0000182' );

    return \%table;
}

sub _do_blast_search {
    my $self = shift;
    my ($fam, $famsize, $fas, $threads, $fam_lengths) = @_;

    my $out = $fas.'.bln';
    my $blastdb = $self->make_blastdb($fas);
    my $blast_report = $self->run_blast({ query   => $fas, 
					  db      => $blastdb, 
					  threads => $threads, 
					  outfile => $out, 
					  evalue  => 1e-10,
					  sort    => 'bitscore' });

    open my $in, '<', $blast_report or die "\n[ERROR]: Could not open file: $blast_report\n";
    return undef unless -e $blast_report && -s $blast_report;

    my (%hits, @pid, @lengths, %uniq);
    while (my $line = <$in>) {
        chomp $line;
        my @f = split /\t/, $line;
        if ($f[3] >= 80) { # hits are at least 80 bp
            push @pid, $f[2];
            push @lengths, $f[3];
            $uniq{$f[0]} = 1;
            $uniq{$f[1]} = 1;
        }
    }
    close $in;
    my @dbfiles = glob "$blastdb*";
    unlink @dbfiles;
    unlink $blast_report;

    # similarity
    my $pid_sum  = sum(@pid);
    my $pid_ct   = @pid;
    my $pid_mean = sprintf("%.2f", $pid_sum/$pid_ct);
    my $count    = keys %uniq;

    # ave HSP length
    my $len_sum  = sum(@lengths);
    my $len_ct   = @lengths;
    my $len_mean = sprintf("%.2f", $len_sum/$len_ct);

    # ave element length
    my $fam_len_sum  = sum(@$fam_lengths);
    my $fam_len_ct   = @$fam_lengths;
    my $fam_len_mean = sprintf("%.2f", $fam_len_sum/$fam_len_ct);

    my $results = join "\t", $fam, $famsize, $count, $pid_mean, $len_mean, $fam_len_mean;

    return $results;
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

    perldoc Tephra::Annotation::Util


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
