package Tephra::Annotation::FilterTandems;

use 5.014;
use Moose;
use File::Spec;
#use File::Find;
use File::Basename;
use File::Path          qw(make_path remove_tree);
use File::Temp          qw(tempfile);
use File::Copy          qw(move);
use Bio::GFF3::LowLevel qw(gff3_parse_feature gff3_format_feature);
use Sort::Naturally;
use Path::Class::File;
#use Data::Dump::Color;
use namespace::autoclean;

with 'Tephra::Role::File',
     'Tephra::Role::Util',
     'Tephra::Role::GFF',
     'Tephra::Role::Run::Blast';

=head1 NAME

Tephra::Role::FilterTandems - Utility methods for filtering LTR/TIR elements that are actually tandemly repeated genes

=head1 VERSION

Version 0.12.5

=cut

our $VERSION = '0.12.5';
$VERSION = eval $VERSION;

has genome => (
    is       => 'ro',
    isa      => 'Path::Class::File',
    required => 1,
    coerce   => 1,
);

has genefile => (
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

has type => (
    is        => 'ro',
    isa       => 'Maybe[Str]',
    predicate => 'has_type',
    required  => 0,
);

has threads => (
    is        => 'ro',
    isa       => 'Int',
    predicate => 'has_threads',
    lazy      => 1,
    default   => 1,
);

#
# methods
#
sub filter_tandem_genes {
    my $self = shift;
    my $gff      = $self->gff->absolute->resolve;
    my $genefile = $self->genefile->absolute->resolve;
    my $genome   = $self->genome->absolute->resolve;
    my $type     = $self->type;

    my $repeat_file = $self->extract_flanking_repeats($gff, $genome, $type);
    my ($filtered_tes, $filtered_stats) = $self->filter_hits_by_percent_overlap($repeat_file, $genefile);
    my $filtered_gff = $self->gff_filter_te_list($gff, $filtered_tes);

    return ($filtered_gff, $filtered_stats);
}

sub extract_flanking_repeats {
    my $self = shift;
    my ($gff, $genome, $type) = @_;

    my $index = $self->index_ref($genome);
    my ($header, $features) = $self->collect_gff_features($gff);

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    
    my ($repeat_files, $repeat_dir, $outfile);

    if ($type eq 'LTR') { 
	$outfile = File::Spec->catfile($path, $name.'_ltrs.fasta');
	open my $outfh, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	
	($repeat_files, $repeat_dir) = $self->extract_ltr_sequences;
	for my $file (@$repeat_files) { 
	    $self->collate($file, $outfh);
	}
	close $outfh;

	remove_tree( $repeat_dir, { safe => 1 } );
	#dd $outfile;
	#dd $repeat_dir and exit;
    }
    elsif ($type eq 'TIR') {
	$outfile = File::Spec->catfile($path, $name.'_tirs.fasta');
	open my $outfh, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";

	($repeat_files, $repeat_dir) = $self->extract_tir_sequences;
	for my $file (@$repeat_files) {
            $self->collate($file, $outfh);
        }
        close $outfh;

	remove_tree( $repeat_dir, { safe => 1 } );
	#dd $repeat_files;
	#dd $repeat_dir and exit;
    }

    return $outfile;
}

sub filter_hits_by_percent_overlap {
    my $self = shift;
    my ($repeat_file, $genefile) = @_;

    my (%removed, %filtered_stats);
    my $repeat_lengths = $self->_store_repeat_lengths($repeat_file);
    my $tect = keys %$repeat_lengths;
    my $rmct = keys %removed;

    my $frac = 0.50; # fraction coverage
    my $pid  = 70;   # percent identity threshold

    my ($rdb, $blast_report, $rdb_is_compressed) = $self->_process_repeat_blast_args($repeat_file, $genefile);
    unless (-s $blast_report) {
	unlink $rdb if $rdb_is_compressed;
	unlink $blast_report;
	unlink $repeat_file;

	my $filtered_tes = $self->_remove_id_store($repeat_lengths, \%removed);

	my $elem_num = $tect / 2;
	$filtered_stats{total}    = $elem_num;
	$filtered_stats{removed}  = $rmct;
	$filtered_stats{filtered} = ($elem_num-$rmct);

	return ($filtered_tes, \%filtered_stats);
    }

    open my $bin, '<', $blast_report or die "\n[ERROR]: Could not open file: $blast_report\n";
    while (my $line = <$bin>) {
	chomp $line;
	my @f = split /\t/, $line;
	if (exists $repeat_lengths->{$f[0]}) {
	    if ($f[3]/$repeat_lengths->{$f[0]} >= $frac) {
		$removed{$f[0]} = sprintf("%.2f",$f[3]/$repeat_lengths->{$f[0]})
		    unless exists $removed{$f[0]};
	    }
	}
	else {
	    warn $f[0]," not found";
	}
    }
    close $bin;

    my $elem_num = $tect / 2;
    $filtered_stats{total}    = $elem_num;
    $filtered_stats{removed}  = $rmct;
    $filtered_stats{filtered} = ($elem_num-$rmct);
    
    my $filtered_tes = $self->_remove_id_store($repeat_lengths, \%removed);
    unlink $rdb if $rdb_is_compressed;
    unlink $blast_report;
    unlink $repeat_file;

    return ($filtered_tes, \%filtered_stats);
}

sub gff_filter_te_list {
    my $self = shift;
    my ($gff, $filtered_tes) = @_;

    #dd $filtered_tes;
    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $tmpiname = 'tephra_transposons_filtered_tandem_genes_XXXX';
    my ($tmp_fh, $tmp_gff) = tempfile( TEMPLATE => $tmpiname, DIR => $path, SUFFIX => '.gff3', UNLINK => 0 );
    my ($header, $features) = $self->collect_gff_features($gff); 
    say $tmp_fh $header;

    my %seen;
    for my $chr_id (nsort keys %$features) {
	for my $rep_region (nsort keys %{$features->{$chr_id}}) {
	    for my $feature (@{$features->{$chr_id}{$rep_region}}) {

		my $teid = $feature->{attributes}{ID}[0];
		next unless defined $teid;

		if ($teid =~ /(?:LTR|TRIM|LARD)_retrotransposon\d+|terminal_inverted_repeat_element\d+|MITE\d+/) {
		    next if exists $seen{$teid};
		    my ($source, $strand) = @{$feature}{qw(source strand)};
		    my $parent_feat = join "||", $chr_id, $source, $rep_region, $strand;
		    if (exists $filtered_tes->{$teid}) {
			$self->_write_features({ parent => $parent_feat, features => $features->{$chr_id}{$rep_region} }, $tmp_fh);
		    }
		    
		    $seen{$teid} = 1;
		}
	    }
	}
    }
    close $tmp_fh;

    #dd $tmp_gff and exit;
    move $tmp_gff, $gff or die "\n[ERROR]: move failed: $!\n";
    
    return $gff->stringify; # this is Path::Class::File object
}

sub _remove_id_store {
    my $self = shift;
    my ($idmap, $removed) = @_;

    my %filtered_tes;
    for my $te (keys %$idmap) {
	my ($teid) = ($te =~ /((?:LTR|TRIM|LARD)_retrotransposon\d+|terminal_inverted_repeat_element\d+|MITE\d+)/);
        unless (exists $removed->{$te}) {
            $filtered_tes{$teid} = 1;
        }
    }

    return \%filtered_tes;
}

sub _process_repeat_blast_args {
    my $self = shift;
    my ($fasta, $rdb) = @_;
    my $threads = $self->threads;

    my $outfile;
    my $rdb_is_compressed = 0;
    if ($rdb =~/\.gz$|\.bz2$/) {
	my $fh = $self->get_fh($rdb);

	$fasta =~ s/\.gz$|\.bz2$//;
        $rdb =~ s/\.gz$|\.bz2$//;
	my ($dbname, $dbpath, $dbsuffix) = fileparse($rdb, qr/\.[^.]*/);
	my ($faname, $fapath, $fasuffix) = fileparse($fasta, qr/\.[^.]*/);
	$outfile = File::Spec->catfile($fapath, $faname.'_'.$dbname.'.bln');
	
	my $tmpiname = $dbname.'_XXXX';
	my ($tmp_fh, $tmp_rdb) = tempfile( TEMPLATE => $tmpiname, DIR => $dbpath, SUFFIX => '', UNLINK => 0 );
	
	while (my $line = <$fh>) {
	    chomp $line;
	    say $tmp_fh $line;
	}
	close $fh;
	close $tmp_fh;

	$rdb = $tmp_rdb;
	$rdb_is_compressed = 1;
    }
    else {
	my ($dbname, $dbpath, $dbsuffix) = fileparse($rdb, qr/\.[^.]*/);
	my ($faname, $fapath, $fasuffix) = fileparse($fasta, qr/\.[^.]*/);

	$outfile = File::Spec->catfile($fapath, $faname.'_'.$dbname.'.bln');
    }

    my $blastdb = $self->make_blastdb($rdb);
    my $blast_report = $self->run_blast({ query   => $fasta, 
                                          db      => $blastdb, 
                                          threads => $threads, 
                                          outfile => $outfile,
                                          sort    => 'bitscore' });
    my @dbfiles = glob "$blastdb*";
    unlink @dbfiles;

    return ($rdb, $blast_report, $rdb_is_compressed);
}

sub _store_repeat_lengths {
    my $self = shift;
    my ($fasta) = @_;

    my %repeats;
    my $kseq = Bio::DB::HTS::Kseq->new($fasta);
    my $iter = $kseq->iterator;
    while (my $seqobj = $iter->next_seq) {
	my $id = $seqobj->name;
	my $seq = $seqobj->seq;
	$repeats{$id} = length($seq);
    }

    return \%repeats;
}

sub _write_features {
    my $self = shift;
    my ($ref, $tmp_fh) = @_;

    my @rep_features = split /\|\|/, $ref->{parent};
    my $parent_feat = 
	join "\t", @rep_features[0,1], 'repeat_region', @rep_features[3,4], '.', $rep_features[5], '.', "ID=$rep_features[2]";
    chomp $parent_feat;
    say $tmp_fh $parent_feat;

    for my $feat (@{$ref->{features}}) {
	my $gff3_str = gff3_format_feature($feat);
	chomp $gff3_str;
	say $tmp_fh $gff3_str;
    }
    say $tmp_fh '###';

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

    perldoc Tephra::Role::FilterTandems


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
