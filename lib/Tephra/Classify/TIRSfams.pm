package Tephra::Classify::TIRSfams;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::SeqIO;
use Bio::Tools::GFF;
use IPC::System::Simple qw(capture);
use Sort::Naturally     qw(nsort);
use List::UtilsBy       qw(nsort_by);
use List::Util          qw(sum max);
use Path::Class::File;
use Try::Tiny;
use Cwd;
use namespace::autoclean;

with 'Tephra::Role::GFF',
     'Tephra::Role::Util';

=head1 NAME

Tephra::Classify::TIRSams - Classify TIR transposons into superfamilies

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has gff => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

#
# methods
#
sub find_tc1_mariner {
    my $self = shift;
    my ($feature, $header) = @_;
    my $fasta = $self->genome;
    my $gff   = $self->gff;
    my $samtools = File::Spec->catfile($ENV{HOME}, '.tephra', 'samtools-1.2', 'samtools');

    my @lengths;
    my $mar_feats;
    my $is_mariner = 0;
    my $has_pdoms  = 0;
    my $pdoms      = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_tc1-mariner.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_tc1-mariner_domain_org.tsv");
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] = $self->_format_attribute($feats[8]);
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=\"?(\w+)\"?/);
		push @all_pdoms, $doms;
	    }
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		if ($tsd_len == 2) {
		    #my $tsd = substr $hash->{$loc}, $feats[3], $tsd_len;
		    my $tmp = $tir.".fasta";
		    my $cmd = "$samtools faidx $fasta $region:$feats[3]-$feats[4] > $tmp";
		    $self->run_cmd($cmd);
		    my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
		    while (my $seqobj = $seqio->next_seq) {
			my $seq = $seqobj->seq;
			#$seq =~ s/.{60}\K/\n/g;
			#say $ofas join "\n", ">".$id, $seq;
			$is_mariner = 1 if $seq =~ /ta/i;
		    }
		    unlink $tmp;
		    
		    #$is_mariner = 1 if $tsd =~ /ta/i;
		}
	    }
	    $mar_feats .= join "\t", @feats, "\n";
	}
	if ($is_mariner) {
	    chomp $mar_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out $mar_feats;
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	}
	undef $mar_feats;
	undef $pdom_org;
	@all_pdoms = ();
	$is_mariner = 0;
	$has_pdoms  = 0;
    }
    close $out;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub find_hat {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $hat_feats;
    my $is_hat = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_hAT.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_hAT_domain_org.tsv");
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";;
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =$self->_format_attribute($feats[8]);
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=\"?(\w+)\"?/);
		push @all_pdoms, $doms;
	    }
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		$is_hat = 1 if $tsd_len == 8;
	    }
	    $hat_feats .= join "\t", @feats, "\n";
	}
	if ($is_hat) {
	    chomp $hat_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out $hat_feats;
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	}
	undef $hat_feats;
	undef $pdom_org;
	@all_pdoms = ();
	$is_hat = 0;
	$has_pdoms = 0;
    }
    close $out;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;    

    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
	
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub find_mutator {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $mut_feats;
    my $is_mutator = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_mutator.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_mutator_domain_org.tsv");
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";;
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =$self->_format_attribute($feats[8]);
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=\"?(\w+)\"?/);
		push @all_pdoms, $doms;
	    }
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		if ($tsd_len >= 8 && $tsd_len <= 11) {
		    $is_mutator = 1;
		}
	    }
	    $mut_feats .= join "\t", @feats, "\n";
	}
	if ($is_mutator) {
	    chomp $mut_feats;
	    say $out $mut_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	}
	undef $mut_feats;
	undef $pdom_org;
	@all_pdoms = ();
	$is_mutator = 0;
	$has_pdoms = 0;
    }
    close $out;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub find_cacta {
    my $self = shift;
    my ($feature, $header) = @_;
    my $fasta = $self->genome;
    my $gff   = $self->gff;
    my $samtools = File::Spec->catfile($ENV{HOME}, '.tephra', 'samtools-1.2', 'samtools');

    my @lengths;
    my $cac_feats;
    my $is_cacta = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_cacta.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_cacta_domain_org.tsv");
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";;
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =$self->_format_attribute($feats[8]);
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=\"?(\w+)\"?/);
		push @all_pdoms, $doms;
	    }
	    if ($feats[2] eq 'target_site_duplication') {
		my $tsd_len = ($feats[4] - $feats[3]) + 1;
		if ($tsd_len >= 2 && $tsd_len <= 3) {
		    #my $tir_elem = substr $hash->{$loc}, $s, $len;
		    my $tmp = $tir.".fasta";
		    my $cmd = "$samtools faidx $fasta $region:$feats[3]-$feats[4] > $tmp";
		    $self->run_cmd($cmd);
		    my $seqio = Bio::SeqIO->new( -file => $tmp, -format => 'fasta' );
		    while (my $seqobj = $seqio->next_seq) {
			my $seq = $seqobj->seq;
			#$seq =~ s/.{60}\K/\n/g;
			#say $ofas join "\n", ">".$id, $seq;
			$is_cacta = 1 if $seq =~ /^cact[ag]/i;
		    }
		    unlink $tmp;
		    
		    #$is_cacta = 1 if $tir_elem =~ /^cact[ag]/i;
		}
	    }
	    $cac_feats .= join "\t", @feats, "\n";
	}
	if ($is_cacta) {
	    chomp $cac_feats;
	    say $out $cac_feats;
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    delete $feature->{$tir};
	    push @lengths, $len;
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	}
	undef $cac_feats;
	undef $pdom_org;
	@all_pdoms = ();
	$is_cacta = 0;
	$has_pdoms = 0;
    }
    close $out;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";
		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub write_unclassified_tirs {
    my $self = shift;
    my ($feature, $header) = @_;
    my $gff = $self->gff;
    
    my @lengths;
    my $unc_feats;
    my $is_unclass = 0;
    my $has_pdoms = 0;
    my $pdoms = 0;
    my %pdom_index;
    my $pdom_org;
    my @all_pdoms;
    
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name."_unclassified.gff3");
    my $domoutfile = File::Spec->catfile($path, $name."_unclassified_domain_org.tsv");
    open my $out, '>>', $outfile or die "\nERROR: Could not open file: $outfile\n";;
    open my $domf, '>>', $domoutfile or die "\nERROR: Could not open file: $domoutfile\n";;
    say $out $header;

    for my $tir (nsort_by { m/repeat_region(\d+)/ and $1 } keys %$feature) {
	my ($rreg, $s, $e) = split /\./, $tir;
	my $len = ($e - $s) + 1;
	my $region = @{$feature->{$tir}}[0];
	my ($loc, $source, $strand) = (split /\|\|/, $region)[0,1,6];
	for my $feat (@{$feature->{$tir}}) {
	    my @feats = split /\|\|/, $feat;
	    $feats[8] =$self->_format_attribute($feats[8]);
	    if ($feats[2] =~ /protein_match/) {
		$has_pdoms = 1;
		my ($doms) = ($feats[8] =~ /name=\"?(\w+)\"?/);
		push @all_pdoms, $doms;
	    }
	    say $out join "\t", $loc, $source, 'repeat_region', $s, $e, '.', '?', '.', "ID=$rreg";
	    say $out join "\t", @feats;
	}
	delete $feature->{$tir};
	push @lengths, $len;
	$pdom_org = join ",", @all_pdoms;
	$pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	$pdoms++ if $has_pdoms;
	$has_pdoms = 0;
    }
    close $out;

    if (%pdom_index) {
	say $domf join "\t", "Strand", "Domain_organizaion", "Domain_count";
	for my $strand (keys %pdom_index) {
	    for my $org (keys %{$pdom_index{$strand}}) {            
		say $domf join "\t", $strand, $org, $pdom_index{$strand}{$org};
	    }
	}
    }
    close $domf;
    unlink $domoutfile unless -s $domoutfile;
    
    my $stat = Statistics::Descriptive::Full->new;
    $stat->add_data(@lengths);
    my $min   = $stat->min;
    my $max   = $stat->max;
    my $mean  = $stat->mean;
    my $count = $stat->count;

    if ($count > 0) {
	say STDERR join "\t", "mariner_count", "min_length", "max_length", "mean_length", "elements_with_protein_matches";		
	say STDERR join "\t", $count, $min, $max, sprintf("%.2f", $mean), $pdoms;
    }
    else {
	unlink $outfile;
    }
}

sub index_ref {
    my $self = shift;
    my $fasta = $self->genome;
    my $samtools  = File::Spec->catfile($ENV{HOME}, '.tephra', 'samtools-1.2', 'samtools');
    my $faidx_cmd = "$samtools faidx $fasta";
    $self->run_cmd($faidx_cmd);
}

sub _format_attribute {
    my $self = shift;
    my ($str) = @_;

    $str =~ s/\s\;\s/\;/g;
    $str =~ s/\s+/=/g;
    $str =~ s/\s+$//;
    $str =~ s/=$//;
    $str =~ s/=\;/;/g;
    $str =~ s/\"//g;
    
    return $str;
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

    perldoc Tephra::Classify::TIRSFams


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
