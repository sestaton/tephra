package Tephra::LTR::Role::Utils;

use 5.014;
use Moose::Role;
use Sort::Naturally;
use File::Spec;
use File::Find;
use File::Basename;
use File::Temp qw(tempfile);
use File::Path qw(remove_tree);
use Bio::DB::HTS::Kseq;
use Carp 'croak';
use namespace::autoclean;
#use Data::Dump::Color;

=head1 NAME

Tephra::LTR::Role::Utils - Common utility methods for working with LTR retrotransposons

=head1 VERSION

Version 0.14.0

=cut

our $VERSION = '0.14.0';
$VERSION = eval $VERSION;

sub get_exemplar_ltrs_for_age {
    my $self = shift;
    my ($dir, $outdir) = @_;

    my (@dirs, @ltrseqs);
    find( sub { push @dirs, $File::Find::name if -d && /_copia\z|_gypsy\z|_unclassified\z/ }, $dir);
    croak "\n[ERROR]: Could not find the expected sub-directories ending in 'copia', 'gypsy' and 'unclassified'. Please ".
        "check input. Exiting.\n" unless @dirs;

    for my $sfdir (@dirs) {
	my ($ltrfile, %ltrfams);
	find( sub { $ltrfile = $File::Find::name if -f and /exemplar_repeats.fasta$/ }, $sfdir);
	unless (defined $ltrfile) {
	    say STDERR "\n[WARNING]: No exemplar LTR file was found in: $sfdir.";
	    say STDERR "This is likely because there were no families identified by the 'classifyltrs' command for this superfamily.";
	    say STDERR "You can try the 'age' command again with the --all flag to process all LTR-RTs.\n";
	    remove_tree( $outdir, { safe => 1 } )
		if $self->clean;
	    exit;
	}

	my $kseq = Bio::DB::HTS::Kseq->new($ltrfile);
	my $iter = $kseq->iterator();
	
	my $re = qr/LTR_retrotransposon\d+|LARD\d+|TRIM\d+/;
	while ( my $seq = $iter->next_seq() ) {
	    my $id  = $seq->name;
	    my $seq = $seq->seq;
	    if ($id =~ /^(?:[35]prime_)?(\w{3}_)?((?:singleton_)?(?:family\d+_))?$re?_\S+_\d+[-_]\d+/) {
		my $code = $1;
		my $family = $2;
		$family =~ s/_//g;
		push @{$ltrfams{$code.$family}}, { id => $id, seq => $seq };
	    }
	}

	for my $family (keys %ltrfams) {
	    my $outfile = File::Spec->catfile($outdir, $family.'_exemplar_ltrseqs.fasta');
	    open my $out, '>', $outfile or die "\n[ERROR]: Could not open file: $outfile\n";
	    for my $pair (@{$ltrfams{$family}}) {
		say $out join "\n", ">".$pair->{id}, $pair->{seq};
	    }
	    close $out;
	    push @ltrseqs, $outfile;
	}
    }

    return \@ltrseqs;
}

sub get_exemplar_ltrs_for_sololtrs {
    my $self = shift;
    my ($solo_obj) = @_;

    my ($dir, $allfams) = @{$solo_obj}{qw(input_dir full_analysis)};
    my ($input, $ltrfile, @ltr_files, @ltrseqs, %ltrfams);
    my $search = qr/[35]prime-ltrs.fasta/;
    my $re = qr/LTR_retrotransposon\d+|LARD\d+|TRIM\d+/;

    if ($allfams) { 
	find( sub { push @ltr_files, $File::Find::name if -f and /$search$/ }, $dir);
    }
    else {
	$search = qr/exemplar_repeats.fasta/;
	find( sub { $ltrfile = $File::Find::name if -f and /$search$/ }, $dir);
	
	if (defined $ltrfile) { 
	    if ($ltrfile =~ /^RL|family\d+/) {
		croak "\n[ERROR]: Expecting a single file of LTR exemplar sequences but it appears this command has ".
		    "been run before. This will cause problems. Please re-run 'classifyltrs' or report this issue. Exiting.\n";
		return;
	    }

	    my $kseq = Bio::DB::HTS::Kseq->new($ltrfile);
	    my $iter = $kseq->iterator();

	    while ( my $seq = $iter->next_seq() ) {
		my $id  = $seq->name;
		my $seq = $seq->seq;
		if ($id =~ /^(?:[35]prime_)?(\w{3}_(?:singleton_)?(?:family\d+_))?($re?)_.*/) {
		    my $family = $1;
		    my $elemid = $2;
		    $family =~ s/^_|_$//;
		    my $name = join "_", $family, $elemid;
		    push @{$ltrfams{$name}}, { id => $id, seq => $seq };
		}
	    }
	}
	else { 
            say STDERR "\n[WARNING]: Exemplar LTR file not found in $dir. This may indicate no families were found. ".
                "Will search all elements instead of exemplars.\n";
            $search = qr/[35]prime-ltrs.fasta/;
            find( sub { push @ltr_files, $File::Find::name if -f and /$search$/ }, $dir);
        }
    }

    if ($allfams || @ltr_files) {
	my (@classified, %name_map);
	find( sub { push @classified, $File::Find::name if -f and /^\w{3}_families.fasta$|^\w{3}_singletons.fasta$/ }, $dir);
	unless (@classified) {
	    croak "\n[ERROR]: Could not find file of LTR singleton or family FASTA files. This indicates an error with the classification. ".
                "Please re-run 'classifyltrs' or report this issue. Exiting.\n";
            return;
	}
	#dd \@classified;

	for my $file (@classified) {
	    my $kseq = Bio::DB::HTS::Kseq->new($file);
            my $iter = $kseq->iterator();
	    
            while ( my $seq = $iter->next_seq() ) {
		my $name;
		my $id = $seq->name;
		my $seq = $seq->seq;
		my ($family, $elemid) = ($id =~ /^(?:[35]prime_)?(\w{3}_(?:singleton_)?(?:family\d+_))?($re?)_.*/); # { 

		if (defined $family && defined $elemid) { 
		    $elemid =~ s/_\d+_\d+$//;
		    $family =~ s/^_|_$//;
		    $name = join "_", $family, $elemid;
		    #$name_map{$name} = $family;
		    $name_map{$elemid} = $family;
		}
		elsif (defined $elemid) {
		    $name_map{$elemid} = 'Unclassified';
		}
		else {
		    $name_map{$id} = 'Unclassified';
		}
            }
	}
	#dd \%name_map;
	
	my $ltr_orient;
	for my $input (nsort @ltr_files) {
	    $ltr_orient = $input =~ /3prime/ ? '3prime' : '5prime';
	    my $kseq = Bio::DB::HTS::Kseq->new($input);
	    my $iter = $kseq->iterator();
	    
	    while ( my $seq = $iter->next_seq() ) {
		my $name;
		my $id = $seq->name;
		my $seq = $seq->seq;
		my ($family, $elemid) = ($id =~ /^(?:[35]prime_)?(\w{3}_(?:singleton_)?(?:family\d+_))?($re?)_.*/);
		if (defined $family && defined $elemid) { 
		    $elemid =~ s/_\d+_\d+$//;
		    $family =~ s/^_|_$//;
		    $name = join "_", $family, $elemid;
		}
		elsif (defined $elemid) { 
		    # No family is defined for this element
		    $name = $elemid;  
		}
		else {
		    $name = $id;
		}

		if (exists $name_map{$name}) {
		    #my $name = join "_", $ltr_orient, $name_map{$elemid}, $id;  
		    my $ltrid = join "_", $ltr_orient, $id;
		    push @{$ltrfams{ $name }}, { id => $ltrid, seq => $seq  }; 
		}
		else {
		    say STDERR "\n[ERROR]: $elemid ($id) not found in family/singleton name map. This is a bug, please report it. Exiting.\n";
		    return;
		}
	    }
	}
    }

    return \%ltrfams;
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

    perldoc Tephra::LTR::Role::Utils


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

1;
