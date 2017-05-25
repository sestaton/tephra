package Tephra::TRIM::TRIMSearch;

use 5.014;
use Moose;
use Path::Class::File;
use File::Spec;
use File::Find;
use File::Basename;
use Cwd                 qw(getcwd abs_path);
use IPC::System::Simple qw(system EXIT_ANY);
use Log::Any            qw($log);
use Try::Tiny;
use namespace::autoclean;

with 'Tephra::Role::Run::GT';

=head1 NAME

Tephra::TRIM::TRIMSearch - Find TRIM retrotransposons in a reference genome

=head1 VERSION

Version 0.07.3

=cut

our $VERSION = '0.07.3';
$VERSION = eval $VERSION;

sub trim_search {
    my $self = shift;
    my ($search_obj) = @_;

    my ($index, $mode) = @{$search_obj}{qw(index mode)};
    my $genome = $self->genome->absolute->resolve;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;

    my (%ltrh_cmd, %ltrd_cmd, @ltrh_opts, @ltrh_args, $ltrh_gff, $ltrg_gff);
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }

    if ($mode eq 'strict') {
	$ltrh_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrharvest99.gff3' );
	$ltrg_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrdigest99.gff3' );

	@ltrh_opts = qw(-seqids -mintsd -maxtsd -minlenltr -maxlenltr 
                       -mindistltr -maxdistltr -motif -similar -vic -index -gff3);
	@ltrh_args = ("yes","4","6","70","500","280","1500","tgca","99","10",$index,$ltrh_gff);
	#@ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    }
    elsif ($mode eq 'relaxed') {
	$ltrh_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrharvest85.gff3' );
	$ltrg_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrdigest85.gff3' );

	@ltrh_opts = qw(-seqids -mintsd -maxtsd -minlenltr -maxlenltr -mindistltr 
                       -maxdistltr -similar -vic -index -gff3);
	@ltrh_args = ("yes","4","6","70","500","280","1500","85","10",$index,$ltrh_gff);
	#@ltrh_cmd{@ltrh_opts} = @ltrh_args;
    }
    else {
	say STDERR "\nERROR: Could not get 'mode' for TRIM search. This is a bug, please report it. Exiting.\n";
	exit(1);
    }

    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff) if -s $ltrh_gff;

    if (defined $gffh_sort && -s $gffh_sort) {
	my @ltrd_opts = qw(-trnas -hmms -seqfile -seqnamelen -matchdescstart -o);
	my @ltrd_args = ($trnadb,$hmmdb,$genome,"50","yes",$ltrg_gff);

	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	$self->clean_indexes($path) 
	    if $self->clean && $mode eq 'relaxed';
	unlink $ltrh_gff;
	unlink $gffh_sort;
    
	return $ltrg_gff;
    }
    else {
	$self->clean_indexes($path) 
	    if $self->clean && $mode eq 'relaxed';
	unlink $ltrh_gff;

	return 0;
    }
}

sub trim_search_strict {
    my $self = shift;
    my ($index) = @_;
    
    my $genome = $self->genome->absolute->resolve;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;

    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrharvest99.gff3' );
    my $ltrg_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrdigest99.gff3' );

    my @ltrh_opts = qw(-seqids -mintsd -maxtsd -minlenltr -maxlenltr 
                       -mindistltr -maxdistltr -motif -similar -vic -index -gff3);
    my @ltrh_args = ("yes","4","6","70","500","280","1500","tgca","99","10",$index,$ltrh_gff);
    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff) if -s $ltrh_gff;

    if (defined $gffh_sort && -s $gffh_sort) {
	my @ltrd_opts = qw(-trnas -hmms -seqfile -seqnamelen -matchdescstart -o);
	my @ltrd_args = ($trnadb,$hmmdb,$genome,"50","yes",$ltrg_gff);
	#my @ltrd_opts = qw(-trnas -hmms -encseq -seqnamelen -matchdescstart -o);
	#my @ltrd_args = ($trnadb,$hmmdb,$index,"50","yes",$ltrg_gff);

	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	$self->clean_indexes($path) if $self->clean;
	#unlink $ltrh_gff;
	#unlink $gffh_sort;
    
	return $ltrg_gff;
    }
    else {
	$self->clean_indexes($path) if $self->clean;
	#unlink $ltrh_gff;
	return 0;
    }
}

sub trim_search_relaxed {
    my $self = shift;
    my ($index) = @_;
    my $genome = $self->genome->absolute->resolve;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;

    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrharvest85.gff3' );
    my $ltrg_gff = File::Spec->catfile( abs_path($path), $name.'_trim_ltrdigest85.gff3' );

    my @ltrh_opts = qw(-seqids -mintsd -maxtsd -minlenltr -maxlenltr -mindistltr 
                       -maxdistltr -similar -vic -index -gff3);
    my @ltrh_args = ("yes","4","6","70","500","280","1500","85","10",$index,$ltrh_gff);
    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff) if -s $ltrh_gff;

    if (defined $gffh_sort && -s $gffh_sort) {
	my @ltrd_opts = qw(-trnas -hmms -seqfile -matchdescstart -seqnamelen -o);
	my @ltrd_args = ($trnadb,$hmmdb,$genome,"yes","50",$ltrg_gff);
	#my @ltrd_opts = qw(-trnas -hmms -encseq -matchdescstart -seqnamelen -o);
	#my @ltrd_args = ($trnadb,$hmmdb,$index,"yes","50",$ltrg_gff);

	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	$self->clean_indexes($path) if $self->clean;
	#unlink $ltrh_gff;
	#unlink $gffh_sort;
	
	return $ltrg_gff;
    }
    else {
	#unlink $ltrh_gff;
	$self->clean_indexes($path) if $self->clean;
	return 0;
    }
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

    perldoc Tephra::TRIM::TRIMSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
