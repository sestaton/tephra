package Tephra::NonLTR::GFFWriter;

use 5.010;
use Moose;
use Bio::SeqIO;
use File::Find;
use File::Spec;
use File::Path qw(make_path);
use File::Basename;
use Sort::Naturally;
use Tephra::Config::Exe;
use namespace::autoclean;

with 'Tephra::Role::Util';

=head1 NAME

Tephra::NonLTR::GFFWriter - Take results from non-LTR search and make an annotated GFF

=head1 VERSION

Version 0.02.1

=cut

our $VERSION = '0.02.1';
$VERSION = eval $VERSION;

has fastadir => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );
has outdir   => ( is => 'ro', isa => 'Maybe[Str]', required => 1 );

sub write_gff {
    my $self = shift;
    my $outdir = $self->outdir;

    my @all_clade = ('CR1', 'I', 'Jockey', 'L1', 'L2', 'R1', 'RandI', 'Rex', 'RTE', 'Tad1', 'R2', 'CRE');
    my $seq_dir   = File::Spec->catdir($outdir, 'info', 'full');
    my @cladedirs = map { File::Spec->catdir($seq_dir, $_) } @all_clade;
    my @nonltrs;

    for my $clade (@cladedirs) {
        my $name = basename($clade);
	my $filename = $name.".dna";
        if (-d $clade) {
	    find( sub { push @nonltrs, $File::Find::name if -f and /$filename/ }, $clade ); 	    
        }
    }

    $self->_fasta_to_gff(\@nonltrs);
    say STDERR "Done with non-LTR search.";
}

sub _fasta_to_gff {
    my $self = shift;
    my $outdir = $self->outdir;
    my ($seqs) = @_;

    my $config = Tephra::Config::Exe->new->get_config_paths;
    my ($samtools) = @{$config}{qw(samtools)};
    my $name = basename($outdir);
    my $outfile = File::Spec->catfile($outdir, $name."_tephra_nonltr.gff3");
    my $fas     = File::Spec->catfile($outdir, $name."_tephra_nonltr.fasta");
    open my $out, '>', $outfile or die "\nERROR: Could not open file: $outfile\n";
    open my $faout, '>', $fas or die "\nERROR: Could not open file: $fas\n";
    my ($lens, $combined) = $self->_get_seq_region;
    $self->_index_ref($samtools, $combined) unless -e $combined.'.fai';

    my %regions;
    for my $file (@$seqs) {
	next if $file =~ /tephra/;
	my ($name, $path, $suffix) = fileparse($file, qr/\.[^.]*/);
	open my $in, '<', $file or die "\nERROR: Could not open file: $file\n";
	while (my $line = <$in>) {
	    chomp $line;
	    if ($line =~ /^>/) {
		my ($id, $start, $end) = ($line =~ /^\>(.*\.?fa(?:s?t?a?))_(\d+)-(\d+)/);
		if (defined $id) {		    
		    $regions{$name}{$id}{$start} = join "||", $start, $end;
		}
		else {
		    warn "\nERROR: Could not parse sequence ID for header: '$line'. ".
			"This is a bug, please report it.\n";;
		}
	    }
	}
	close $in;
    }

    my @ids = map { nsort keys %{$regions{$_}} } keys %regions;
    say $out '##gff-version 3';

    for my $seqid (@ids) {
	my ($name, $path, $suffix) = fileparse($seqid, qr/\.[^.]*/);
	if (exists $lens->{$name}) {
	    say $out join q{ }, '##sequence-region', $name, '1', $lens->{$name};
	}
	else {
	    say "\nERROR: Could not find $name in map.\n";
	}
    }

    ##TODO: How to get the strand correct?
    my $ct = 0;
    for my $clade (keys %regions) {
	for my $seqid (nsort keys %{$regions{$clade}}) {
	    my ($name, $path, $suffix) = fileparse($seqid, qr/\.[^.]*/);
	    for my $start (sort { $a <=> $b } keys %{$regions{$clade}{$seqid}}) {
	        my ($start, $end) = split /\|\|/, $regions{$clade}{$seqid}{$start};
		my $seqname = $seqid;
		$seqname =~ s/\.fa.*//;
		my $elem = "non_LTR_retrotransposon$ct";
                my $tmp = $elem.".fasta";
                my $id  = $elem."_".$seqname."_".$start."_".$end;
                my $cmd = "$samtools faidx $combined $seqname:$start-$end > $tmp";
                $self->run_cmd($cmd);
		my $seqio = Bio::SeqIO->new(-file => $tmp, -format => 'fasta');
		while (my $seqobj = $seqio->next_seq) {
		    my $seq = $seqobj->seq;
		    $seq =~ s/.{60}\K/\n/g;
		    say $faout join "\n", ">".$id, $seq;
		}
		unlink $tmp;
		
		say $out join "\t", $name, 'Tephra', 'non_LTR_retrotransposon', $start, $end, '.', '?', '.', 
		    "ID=non_LTR_retrotransposon$ct;Name=$clade;Ontology_term=SO:0000189"; 
		$ct++;
	    }
	} 
    }
    close $faout;
    unlink $combined;
}

sub _get_seq_region {
    my $self = shift;
    my $fasdir = $self->fastadir;

    my (@seqs, %lens);
    find( sub { push @seqs, $File::Find::name if -f and /\.fa.*$/ }, $fasdir );
    my $combined = File::Spec->catfile($fasdir, 'tephra_all_genome_seqs.fas');
    open my $out, '>', $combined or die "\nERROR: Could not open file: $combined\n";

    for my $seq (nsort @seqs) {
	next if $seq =~ /tephra/;
	$self->_collate($seq, $out);
	my $seqio = Bio::SeqIO->new(-file => $seq, -format => 'fasta');
	while (my $seqobj = $seqio->next_seq) {
	    my $id  = $seqobj->id;
	    my $len = $seqobj->length;
	    $lens{$id} = $len;
	}
    }

    return (\%lens, $combined);
}

sub _collate {
    my $self = shift;
    my ($file_in, $fh_out) = @_;
    my $lines = do { 
        local $/ = undef; 
        open my $fh_in, '<', $file_in or die "\nERROR: Could not open file: $file_in\n";
        <$fh_in>;
    };
    print $fh_out $lines;
}

sub _index_ref {
    my $self = shift;
    my ($samtools, $fasta) = @_;
    my $faidx_cmd = "$samtools faidx $fasta";
    $self->run_cmd($faidx_cmd);
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

    perldoc Tephra::NonLTR::GFFWriter


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
