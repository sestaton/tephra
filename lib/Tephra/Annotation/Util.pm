package Tephra::Annotation::Util;

use 5.014;
use Moose;
use List::Util qw(sum);
use Sort::Naturally;
use Bio::DB::HTS::Kseq;
use namespace::autoclean;

with 'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Annotation::Util - Utility methods for working with transposon annotations

=head1 VERSION

Version 0.09.0

=cut

our $VERSION = '0.09.0';
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
    #my $usage = "$0 ltr-fams_combined.fasta similarity.tsv\n";
    my ($fasta, $outfile, $threads) = @_;
    
    open my $outfh, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    my $kseq = Bio::DB::HTS::Kseq->new($fasta);
    my $iter = $kseq->iterator;
    
    my ($famname, %families, %lengths);
    while (my $obj = $iter->next_seq) {
	my $id = $obj->name;
	my $seq = $obj->seq;
	if ($id =~ /^(\w{3}_(?:singleton_)?family\d+)_/) {
	    $famname = $1;
	    push @{$families{$famname}}, { seq => $seq, id => $id };
	    #$lengths{$id} = length($seq);
	}
	else {
	    say STDERR "\nWARNING: '$id' does not match pattern.\n";
	}
    }
    #dd \%families and exit;

    for my $fam (reverse sort { @{$families{$a}} <=> @{$families{$b}} } keys %families) {
	my $famsize = @{$families{$fam}};
	#say join "\t", $fam, $famsize;
	if (@{$families{$fam}} > 1) {
	    my $outfile = $fam.'.fasta';
	    open my $out, '>', $outfile;
	    for my $seqobj (@{$families{$fam}}) {
		say $out join "\n", ">".$seqobj->{id}, $seqobj->{seq};
	    }
	    close $out;
	    $self->_do_blast_search($fam, $famsize, $outfile, $outfh, $threads);
	    unlink $outfile;
	}
    }
    #dd \%families;

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
	'RII' => 'I',
        'RIJ' => 'Jockey',
        'RIL' => 'L1',
        'RIR' => 'R2',
        'RIT' => 'RTE',
	'RIC' => 'CR1',
	'RIS' => 'R1',
	'RIA' => 'RandI',
	'RIE' => 'Rex',
	'RID' => 'Tad1',
        'RIX' => 'Unknown_LINE',
        'RLB' => 'Bel/Pao',
        'RLG' => 'Gypsy',
        'RLC' => 'Copia',
        'RLE' => 'ERV',
        'RLR' => 'Retrovirus',
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
	say STDERR "\n[WARNING]: No 3-letter code could be found for: $sfamily_code\n"
	    if $self->debug;
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
        'RII' => 'I',
        'RIJ' => 'Jockey',
        'RIL' => 'L1',
        'RIR' => 'R2',
        'RIT' => 'RTE',
	'RIC' => 'CR1',
	'RIS' => 'R1',
	'RIA' => 'RandI',
	'RIE' => 'Rex',
	'RID' => 'Tad1',
        'RIX' => 'Unknown_LINE',
        'RLB' => 'Bel/Pao',
        'RLG' => 'Gypsy',
        'RLC' => 'Copia',
        'RLE' => 'ERV',
        'RLR' => 'Retrovirus',
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
        'DTX' => 'Unknown_TIR',
        'DXX' => 'Unknown_DNA_transposon',
        'DHH' => 'Helitron',
        'DHX' => 'Unknown_Helitron',
        'DMM' => 'Maverick',
        'DMX' => 'Unknown_Maverick',
	'RST' => 'SINE2/tRNA' );
    
    my %name_table = reverse %sfcode_table;

    if (exists $name_table{$name}) {
	return $name_table{$name};
    }
    else {
	say STDERR "\n[WARNING]: No 3-letter code could be found for: $name\n"
	    if $self->debug;
	return 0;
    }
}

sub _do_blast_search {
    my $self = shift;
    my ($fam, $famsize, $fas, $outfh, $threads) = @_;

    my $out = $fas.'.bln';
    my $blastdb = $self->make_blastdb($fas);
    #system("makeblastdb -in $fas -dbtype nucl 2>&1 > /dev/null") 
        #== 0 or die $!;
    #system('blastn', '-query', $fas, '-db', $fas, '-outfmt', '6', '-num_threads', 12, '-out', $out, '-evalue', 1e-10) 
        #== 0 or die $!;
    my $blast_report = $self->run_blast({ query   => $fas, 
					  db      => $blastdb, 
					  threads => $threads, 
					  outfile => $out, 
					  evalue  => 1e-10,
					  sort    => 'bitscore' });

    open my $in, '<', $blast_report or die "\nERROR: Could not open file: $blast_report\n";

    my (%hits, @pid, @lengths, %uniq);
    while (my $line = <$in>) {
        chomp $line;
        my @f = split /\t/, $line;
        if ($f[3] >= 80) { # hits are at least 80 bp
            #say $outfh join "\t", $fam, @f[0..3];
            #push @{$hits{$fam}}, { stats => join "||", @f[0..3], pid => $f[3] };
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

    #my $stat = Statistics::Descriptive::Full->new;
    #$stat->add_data(@pid);
    my $pid_sum = sum(@pid);
    my $pid_ct  = @pid;
    my $pid_mean = sprintf("%.2f", $pid_sum/$pid_ct);
    my $count = keys %uniq;
    #my $mean = $stat->mean;

    #my $lstat = Statistics::Descriptive::Full->new;
    #$lstat->add_data(@lengths);
    #my $lmean = $lstat->mean;
    my $len_sum = sum(@lengths);
    my $len_ct  = @lengths;
    my $len_mean = sprintf("%.2f", $len_sum/$len_ct);

    say $outfh join "\t", 'family', 'family_size', 'family_count_analyzed', 'similarity_mean', 'length_mean';
    say $outfh join "\t", $fam, $famsize, $count, $pid_mean, $len_mean;

    return;
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

    perldoc Tephra::Annotation::Util


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
