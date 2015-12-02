package Tephra::LTR::LTRSearch;

use 5.010;
use Moose;
use Cwd;
use File::Spec;
use File::Find;
use File::Basename;
use IPC::System::Simple qw(system EXIT_ANY);
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;
use namespace::autoclean;
use Data::Dump;

with 'Tephra::Role::Run::GT';

=head1 NAME

Tephra::LTR::LTRSearch - Find LTR retrotransposons in a reference genome

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';
$VERSION = eval $VERSION;

has mintsd => ( is => 'ro', isa => 'Int', required => 0, default => 4 );
has maxtsd => ( is => 'ro', isa => 'Int', required => 0, default => 6 );
has minlenltr => ( is => 'ro', isa => 'Int', required => 0, default => 100 );
has maxlenltr => ( is => 'ro', isa => 'Int', required => 0, default => 6000 );
has mindistltr => ( is => 'ro', isa => 'Int', required => 0, default => 1500 );
has maxdistltr => ( is => 'ro', isa => 'Int', required => 0, default => 25000 );

sub ltr_search_strict {
    my $self = shift;
    my ($index) = @_;
    
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;

    my $mintsd = $self->mintsd;
    my $maxtsd = $self->maxtsd;
    my $minlenltr = $self->minlenltr;
    my $maxlenltr = $self->maxlenltr;
    my $mindistltr = $self->mindistltr;
    my $maxdistltr = $self->maxdistltr;

    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_ltrharvest99_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_ltrharvest99_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_ltrharvest99.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_ltrdigest99.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_ltrdigest99");

    my @ltrh_opts = qw(-longoutput -seqids -tabout -mintsd -maxtsd -minlenltr -maxlenltr -mindistltr 
                       -maxdistltr -motif -similar -vic -index -out -outinner -gff3);

    my @ltrh_args = ("no","yes","no",$mintsd,$maxtsd,$minlenltr,$maxlenltr,$mindistltr,$maxdistltr,"tgca","99",
		     "10",$index,$ltrh_out,$ltrh_out_inner,$ltrh_gff);
    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    if (-s $ltrh_gff) {
	my $gffh_sort = $self->sort_gff($ltrh_gff);

        #    -pdomevalcutoff    global E-value cutoff for pHMM search
        #	default 1E-6
        #	-pdomcutoff

	my @ltrd_opts = qw(-trnas -hmms -aliout -aaout -seqfile -matchdescstart -seqnamelen -o -outfileprefix);
        #-pdomevalcutoff -pdomcutoff);
	my @ltrd_args = ($trnadb,$hmmdb,"no","no",$genome,"yes","50",$ltrg_gff,$ltrg_out); #,'1E-10','TC');
	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	unlink $ltrh_gff;
	unlink $gffh_sort;
    
	return $ltrg_gff;
    }
    else {
	unlink $ltrh_gff;
    }
}

sub ltr_search_relaxed {
    my $self = shift;
    my ($index) = @_;
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;
   
    my $mintsd = $self->mintsd;
    my $maxtsd = $self->maxtsd;
    my $minlenltr = $self->minlenltr;
    my $maxlenltr = $self->maxlenltr;
    my $mindistltr = $self->mindistltr;
    my $maxdistltr = $self->maxdistltr;
    
    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_ltrharvest85_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_ltrharvest85_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_ltrharvest85.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_ltrdigest85.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_ltrdigest85");

    my @ltrh_opts = qw(-longoutput -seqids -tabout -mintsd -maxtsd -minlenltr -maxlenltr 
                       -mindistltr -maxdistltr -similar -vic -index -out -outinner -gff3);

    my @ltrh_args = ("no","yes","no",$mintsd,$maxtsd,$minlenltr,$maxlenltr,$mindistltr,$maxdistltr,"85","10",
		     $index,$ltrh_out,$ltrh_out_inner,$ltrh_gff);
    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff);

    my @ltrd_opts = qw(-trnas -hmms -aliout -aaout -seqfile -matchdescstart -seqnamelen -o -outfileprefix);
    my @ltrd_args = ($trnadb,$hmmdb,"no","no",$genome,"yes","50",$ltrg_gff,$ltrg_out);
    @ltrd_cmd{@ltrd_opts} = @ltrd_args;
    
    my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
    $self->clean_index if $self->clean;
    unlink $ltrh_gff;
    unlink $gffh_sort;

    return $ltrg_gff;
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

    perldoc Tephra::LTR::LTRSearch


=head1 LICENSE AND COPYRIGHT

Copyright (C) 2015- S. Evan Staton

This program is distributed under the MIT (X11) License, which should be distributed with the package. 
If not, it can be found here: L<http://www.opensource.org/licenses/mit-license.php>

=cut

__PACKAGE__->meta->make_immutable;

1;
