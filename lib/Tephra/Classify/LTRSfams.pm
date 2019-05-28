package Tephra::Classify::LTRSfams;

use 5.014;
use Moose;
use MooseX::Types::Path::Class;
use Statistics::Descriptive;
use File::Spec;
use File::Find;
use File::Basename;
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use IPC::System::Simple qw(capture);
use List::UtilsBy       qw(nsort_by);
use Cwd                 qw(getcwd abs_path);
use Path::Class         qw(file);
use Sort::Naturally     qw(nsort);
use Try::Tiny;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::Logger',
     'Tephra::Role::GFF',
     'Tephra::Role::Util',
     'Tephra::Role::Run::Blast',
     'Tephra::Classify::Role::LogResults';

=head1 NAME

Tephra::Classify::LTRSFams - Classify LTR retrotransposons into superfamilies

=head1 VERSION

Version 0.12.4

=cut

our $VERSION = '0.12.4';
$VERSION = eval $VERSION;

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has repeatdb => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
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
      required => 1,
      coerce   => 1,
);

has blast_hit_length => (
    is      => 'ro',
    isa     => 'Int',
    default => 80,
);

has blast_hit_pid => (
    is      => 'ro',
    isa     => 'Int',
    default => 80,
);

#
# methods
#
sub find_gypsy_copia {
    my $self = shift;
    my ($features) = @_;

    my ($is_gypsy, $is_copia) = (0, 0);
    my (%gypsy, %copia, %pdom_index, $pdom_org, @all_pdoms, @gyp_cop_doms);

    my @cop_exp = qw(gag asp rve UBN2 UBN2_2 UBN2_3 RVT_2 RNase_H);
    my @gyp_exp = qw(gag asp RVT_1 RNash_H rve Chromo);
    
    my $strand;
    for my $seq_id (keys %$features) {
	for my $rep_region (keys %{$features->{$seq_id}}) {
	    for my $ltr_feature (@{$features->{$seq_id}{$rep_region}}) {
		$strand = $ltr_feature->{strand};
		if ($ltr_feature->{type} eq 'protein_match') {
		    my $pdom_name = $ltr_feature->{attributes}{name}[0];
		    push @all_pdoms, $pdom_name;
		}
	    }
	    if (@all_pdoms) {
		$pdom_org = join ",", @all_pdoms;
		@all_pdoms = reverse @all_pdoms if $strand eq '-';
		# gypsy -> gag,ap,int,rt,rh
		# copia -> gag,ap,rt,rh,int
		#
		# model names: gag|retrotransposon_gag
		# NAME  UBN2_2
		# DESC  gag-polypeptide of LTR copia-type
		# NAME  UBN2_3
		
		# DESC  gag-polypeptide of LTR copia-type
		# NAME  UBN2
		# DESC  gag-polypeptide of LTR copia-type
		# REF   The nucleotide sequence of Drosophila melanogaster 
		#       copia-specific 2.1-kb mRNA.Nucleic Acids Res. 1989 Mar 11; 17(5):2134

		# NAME  Retrotrans_gag
		# DESC  Retrotransposon gag protein
		
		# NAME  gag-asp_proteas
		# DESC  gag-polyprotein putative aspartyl protease
		
		# NAME  gag_pre-integrs
		# DESC  GAG-pre-integrase domain
		# NAME  rve
		# DESC  Integrase core domain
		
		# NAME  RVT_1
		# DESC  Reverse transcriptase (RNA-dependent DNA polymerase)
		# NAME  RVT_2
		# DESC  Reverse transcriptase (RNA-dependent DNA polymerase)
		# NAME  RVT_3
		# DESC  Reverse transcriptase-like
		
		# NAME  RNase_H
		# DESC  RNase H
		
		# NAME  Chromo
		# DESC  Chromo (CHRromatin Organisation MOdifier) domain
		#say STDERR join q{ }, "strand: $strand", $pdom_org;
		my ($gyp_org, $cop_org);
		my ($gyp_dom_ct, $cop_dom_ct) = (0, 0);
		if (grep { /rvt_2|ubn2/i && ! /rvt_1|chromo/i } @all_pdoms) {
		    for my $d (@cop_exp) {
			for my $p (@all_pdoms) {
			    $cop_dom_ct++ if $p =~ /$d/i;
			}
		    }
		}
		
		if (grep { /rvt_1|chromo/i && ! /rvt_2|ubn2/i } @all_pdoms) {
		    for my $d (@gyp_exp) {
			for my $p (@all_pdoms) {
			    $gyp_dom_ct++ if $p =~ /$d/i;
			}
		    }
		}
		
		if ($gyp_dom_ct >= 1) {
		    $gypsy{$seq_id}{$rep_region} = $features->{$seq_id}{$rep_region};
		    push @gyp_cop_doms, join "||", $seq_id, $rep_region;
		}
		elsif ($cop_dom_ct >= 1) {
		    $copia{$seq_id}{$rep_region} = $features->{$seq_id}{$rep_region};
		    push @gyp_cop_doms, join "||", $seq_id, $rep_region;
		}
		
		$gyp_dom_ct = 0;
		$cop_dom_ct = 0;
		$is_gypsy   = 0;
		$is_copia   = 0;
		
		undef $pdom_org;
		@all_pdoms = ();
	    }
	}    
    }

    $self->remove_repeat_region_features($features, \@gyp_cop_doms);

    #dd \%gypsy and exit;
    return (\%gypsy, \%copia);
}

sub find_unclassified {
    my $self = shift;
    my ($features) = @_;
    my $gff   = $self->gff->absolute->resolve;
    my $fasta = $self->genome->absolute->resolve;
    my $index = $self->index_ref($fasta);

    my %ltr_rregion_map;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfast = File::Spec->catfile($path, $name.'_unclassified.fasta');

    open my $ofas, '>>', $outfast or die "\n[ERROR]: Could not open file: $outfast\n";

    for my $seq_id (keys %$features) {
	for my $rep_region (keys %{$features->{$seq_id}}) {
	    for my $ltr_feature (@{$features->{$seq_id}{$rep_region}}) {
		if ($ltr_feature->{type} =~ /(?:LTR|TRIM|LARD)_retrotransposon/) {
		    my ($seq_id, $start, $end) = @{$ltr_feature}{qw(seq_id start end)};
		    my $elem = $ltr_feature->{attributes}{ID}[0];
		    my $id = join "_", $elem, $seq_id, $start, $end;
		    $ltr_rregion_map{$id} = join "||", $seq_id, $rep_region;
		    
		    my ($seq, $length) = $self->get_full_seq($index, $seq_id, $start, $end);
		    $seq =~ s/.{60}\K/\n/g;
		    say $ofas join "\n", ">".$id, $seq;
		}
	    }
	}
    }
    close $ofas;

    return ($outfast, \%ltr_rregion_map);
}

sub search_unclassified {
    my $self = shift;
    my ($unc_fas) = @_;
    my $repeatdb = $self->repeatdb->absolute->resolve;
    my $threads  = $self->threads;
    my $blastdb = $self->make_blastdb($repeatdb);

    my ($bname, $bpath, $bsuffix) = fileparse($repeatdb, qr/\.[^.]*/);
    my ($fname, $fpath, $fsuffix) = fileparse($unc_fas, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($fpath, $fname.'_'.$bname.'.bln');
    my $blast_report = $self->run_blast({ query   => $unc_fas, 
					  db      => $blastdb, 
					  threads => $threads, 
					  outfile => $outfile, 
					  sort    => 'bitscore' });

    my @dbfiles = glob "$blastdb*";
    unlink @dbfiles;

    return $outfile;
}

sub annotate_unclassified {
    my $self = shift;
    my ($blast_out, $gypsy, $copia, $features, $ltr_rregion_map) = @_;
    my $hit_length = $self->blast_hit_length;
    my $hit_pid    = $self->blast_hit_pid;

    my $family_map = $self->_map_repeat_types();
    open my $in, '<', $blast_out or die "\n[ERROR]: Could not open file: $blast_out\n";
    my (%gypsy_re, %copia_re, %seen, @gyp_cop_regions);
    #dd $ltr_rregion_map and exit;

    while (my $line = <$in>) {
	chomp $line;
	my @f = split /\t/, $line;
	next if exists $seen{$f[0]}; # this is taking the best hit only, by bitscore

	if ($f[2] >= $hit_pid && $f[3] >= $hit_length) {
	    if (exists $family_map->{$f[1]}) {
		my ($seq_id, $rregion, $start, $end) = split /\|\|/, $ltr_rregion_map->{$f[0]};
		my $key = join "||", $rregion , $start, $end;
		my $sf = $family_map->{$f[1]};
		if ($sf =~ /^rlg|gypsy/i) {
		    $gypsy->{ $seq_id }{ $key } = $features->{ $seq_id }{ $key };
		    push @gyp_cop_regions, join "||", $seq_id, $key;
		    $seen{$f[0]} = 1;
		}
		elsif ($sf =~ /^rlc|copia/i) {
		    $copia->{ $seq_id }{ $key } = $features->{ $seq_id }{ $key };
		    push @gyp_cop_regions, join "||", $seq_id, $key;
		    $seen{$f[0]} = 1;
		}
	    }
	}
    }
    close $in;
    unlink $blast_out;

    $self->remove_repeat_region_features($features, \@gyp_cop_regions);

    return;
}

sub write_gypsy {
    my $self = shift;
    my ($gypsy, $header, $log) = @_;
    my $gff = $self->gff->absolute->resolve;

    my (@lengths, @all_pdoms);
    my ($gyp_feats, %pdom_index, $pdom_org);
    my ($has_pdoms, $pdoms) = (0, 0);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile    = File::Spec->catfile($path, $name.'_gypsy.gff3');
    my $domoutfile = File::Spec->catfile($path, $name.'_gypsy_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    say $out $header;
    
    my $gyp_ct = (keys %$gypsy);
    #dd $gypsy and exit;

    my ($seq_id, $source, $start, $end, $strand);
    for my $chr_id (nsort keys %$gypsy) {
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$gypsy->{$chr_id}}) {
	    my ($rreg, $s, $e) = split /\|\|/, $rep_region;
	    for my $ltr_feature (@{$gypsy->{$chr_id}{$rep_region}}) {
		if ($ltr_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    my $pdom_name = $ltr_feature->{attributes}{name}[0];
		    push @all_pdoms, $pdom_name;
		}
		if ($ltr_feature->{type} =~ /(?:LTR|TRIM|LARD)_retrotransposon/) {
		    ($seq_id, $source, $start, $end, $strand) 
			= @{$ltr_feature}{qw(seq_id source start end strand)};
		    $strand //= '?';
		    my $ltrlen = $end - $start + 1;
		    push @lengths, $ltrlen;
		}
		my $gff3_str = gff3_format_feature($ltr_feature);
		$gyp_feats .= $gff3_str;
	    }

	    chomp $gyp_feats;
	    say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg";
	    say $out $gyp_feats;
	    
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	    undef $gyp_feats;
	    undef $pdom_org;
	    @all_pdoms = ();
	    $has_pdoms = 0;
	}
    }
    close $out;
    
    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile }) 
	if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) {
	my $count = $self->log_basic_element_stats({ lengths => \@lengths, type => 'Gypsy', log => $log, pdom_ct => $pdoms });

	return $outfile;
    }
    else {
	return undef;
    }
}

sub write_copia {
    my $self = shift;
    my ($copia, $header, $log) = @_;
    my $gff = $self->gff->absolute->resolve;
    
    my (@lengths, $cop_feats);
    my ($pdom_org, @all_pdoms, %pdom_index);
    my ($has_pdoms, $pdoms) = (0, 0);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_copia.gff3');
    my $domoutfile = File::Spec->catfile($path, $name.'_copia_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    say $out $header;

    my ($seq_id, $source, $start, $end, $strand);
    for my $chr_id (nsort keys %$copia) {
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$copia->{$chr_id}}) {
	    my ($rreg, $s, $e) = split /\|\|/, $rep_region;
	    for my $ltr_feature (@{$copia->{$chr_id}{$rep_region}}) {
		if ($ltr_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    my $pdom_name = $ltr_feature->{attributes}{name}[0];
		    push @all_pdoms, $pdom_name;
		}
		if ($ltr_feature->{type} =~ /(?:LTR|TRIM|LARD)_retrotransposon/) {
		    ($seq_id, $source, $start, $end, $strand) 
			= @{$ltr_feature}{qw(seq_id source start end strand)};
		    my $ltrlen = $end - $start + 1;
		    push @lengths, $ltrlen;
		}
		my $gff3_str = gff3_format_feature($ltr_feature);
		$cop_feats .= $gff3_str;
	    }
	    chomp $cop_feats;
	    say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg";
	    say $out $cop_feats;
	    
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	    undef $cop_feats;
	    undef $pdom_org;
	    @all_pdoms = ();
	    $has_pdoms = 0;
	}
    }
    close $out;

    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
        if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;

    if (@lengths) {
	my $count = $self->log_basic_element_stats({ lengths => \@lengths, type => 'Copia', log => $log, pdom_ct => $pdoms });

	return $outfile;
    }
    else {
	return undef;
    }
}

sub write_unclassified {
    my $self = shift;
    my ($features, $header, $log) = @_;
    my $gff = $self->gff->absolute->resolve;

    my (%pdom_index, %lard_index);
    my (@all_pdoms, $pdom_org);
    my (@lard_lengths, @unc_lengths);
    my ($has_pdoms, $pdoms, $is_lard, $lard_ct) = (0, 0, 0, 1);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $outfile = File::Spec->catfile($path, $name.'_unclassified.gff3');
    my $domoutfile = File::Spec->catfile($path, $name.'_unclassified_domain_org.tsv');
    open my $out, '>>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
    say $out $header;

    my ($seq_id, $source, $start, $end, $strand, $ltrlen, $unc_feats);
    for my $chr_id (nsort keys %$features) {
	for my $rep_region (nsort_by { m/repeat_region\d+\.?\d+?\|\|(\d+)\|\|\d+/ and $1 } keys %{$features->{$chr_id}}) {
	    my ($rreg, $s, $e) = split /\|\|/, $rep_region;
	    for my $ltr_feature (@{$features->{$chr_id}{$rep_region}}) {
		if ($ltr_feature->{type} eq 'protein_match') {
		    $has_pdoms = 1;
		    my $pdom_name = $ltr_feature->{attributes}{name}[0];
		    push @all_pdoms, $pdom_name;
		}
		if ($ltr_feature->{type} =~ /(?:LTR|TRIM)_retrotransposon/) {
		    ($seq_id, $source, $start, $end, $strand) 
			= @{$ltr_feature}{qw(seq_id source start end strand)};
		    $strand //= '?';
		    $ltrlen = $end - $start + 1;
		}
		my $gff3_str = gff3_format_feature($ltr_feature);
		$unc_feats .= $gff3_str;
	    }

	    chomp $unc_feats;
	    say $out join "\t", $seq_id, $source, 'repeat_region', $s, $e, '.', $strand, '.', "ID=$rreg";

	    if ($has_pdoms) { 
		say $out $unc_feats;
		push @unc_lengths, $ltrlen;
	    }
	    else {
		my $unc_lard_feats = $self->_format_lard_features($unc_feats, $lard_ct);
		if ($unc_lard_feats->{is_lard}) {
		    $lard_ct++;
		    say $out $unc_lard_feats->{unc_feats};
		    push @lard_lengths, $ltrlen;
		    $lard_index{ $unc_lard_feats->{old_id} } = $unc_lard_feats->{new_id};
		}
		else {
		    say $out $unc_feats;
		    push @unc_lengths, $ltrlen;
		}
	    }
	    
	    $pdom_org = join ",", @all_pdoms;
	    $pdom_index{$strand}{$pdom_org}++ if $pdom_org;
	    $pdoms++ if $has_pdoms;
	    @all_pdoms = ();
	    
	    $has_pdoms = 0;
	    undef $unc_feats;
	    undef $pdom_org;
	}
    }
    close $out;
	
    $self->write_superfam_pdom_organization({ pdom_index => \%pdom_index, outfile => $domoutfile })
	if %pdom_index;
    unlink $domoutfile unless -s $domoutfile;
    
    if (@unc_lengths || @lard_lengths) {
	my $unc_count = 
	    $self->log_basic_element_stats({ lengths => \@unc_lengths, type => 'unclassified LTR-RT', log => $log, pdom_ct => $pdoms });
	my $lard_count =
	    $self->log_basic_element_stats({ lengths => \@lard_lengths, type => 'LARD', log => $log, pdom_ct => 0 });
	
	return ($outfile, \%lard_index);
    }
    else { 
	return undef;
    }
}

sub _format_lard_features {
    my $self = shift;
    my ($unc_feats, $lard_ct) = @_;

    my ($new_feats, $old_id);
    my $is_lard = 0;
    my $lard_type = 'LARD_retrotransposon';
    my $new_id = $lard_type.$lard_ct;

    for my $feat (split /^/, $unc_feats) {
        chomp $feat;
        my @feats = split /\t/, $feat;
	if ($feats[2] eq 'LTR_retrotransposon') {
            $feats[2] = $lard_type;
            $is_lard = 1;
        }
	if ($feats[8] =~ /(LTR_retrotransposon\d+)/) {
	    $old_id = $1;
	    $feats[8] =~ s/$old_id/$new_id/g;
        }
        $new_feats .= join "\t", @feats, "\n";
    }
    chomp $new_feats;

    return ({ unc_feats => $new_feats, is_lard => $is_lard, old_id => $old_id, new_id => $new_id });
}

sub _map_repeat_types {
    my $self = shift;
    my $repeatdb = $self->repeatdb->absolute->resolve;
    my %family_map;

    open my $in, '<', $repeatdb or die "\n[ERROR]: Could not open file: $repeatdb\n";

    while (my $line = <$in>) {
        chomp $line;
        if ($line =~ /^>/) {
            $line =~ s/>//;
            my ($f, $sf, $source)  = split /\t/, $line;
            next unless defined $sf && defined $f;
            if ($sf =~ /(\s+)/) {
                $sf =~ s/$1/\_/;
            }
            $f =~ s/\s/\_/;
            $family_map{$f} = $sf;
        }
    }
    close $in;

    return \%family_map;
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

    perldoc Tephra::Classify::LTRSFams


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut 

__PACKAGE__->meta->make_immutable;

1;
