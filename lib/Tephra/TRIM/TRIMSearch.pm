package Tephra::TRIM::TRIMSearch;

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

with 'Tephra::Role::GT';

sub trim_search_strict {
    my $self = shift;
    my ($index) = @_;
    
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    #my $gt = $self->get_gt_exec;
    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_trim_ltrharvest99_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_trim_ltrharvest99_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_trim_ltrharvest99.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_trim_ltrdigest99.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_trim_ltrdigest99");

    my @ltrh_opts = qw(-seqids -longoutput -mintsd -maxtsd -minlenltr -maxlenltr 
                       -mindistltr -maxdistltr -motif -similar -vic -index -out -outinner -gff3);

    my @ltrh_args = ("yes","yes","4","6","70","500","280","1500","tgca","99","10",
		     $index,$ltrh_out,$ltrh_out_inner,$ltrh_gff);

    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff) if -s $ltrh_gff;

    if (defined $gffh_sort && -s $gffh_sort) {
	my @ltrd_opts = qw(-trnas -hmms -aliout -aaout -seqfile -seqnamelen -matchdescstart -o -outfileprefix);
	my @ltrd_args = ($trnadb,$hmmdb,"no","no",$genome,"50","yes",$ltrg_gff,$ltrg_out);
	
	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	#$self->clean_index if $self->clean;
	unlink $ltrh_gff;
	unlink $gffh_sort;
    
	return $ltrg_gff;
    }
    else {
	return 0;
    }
}

sub trim_search_relaxed {
    my $self = shift;
    my ($index) = @_;
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;
    my (%suf_args, %ltrh_cmd, %ltrd_cmd);
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_trim_ltrharvest85_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_trim_ltrharvest85_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_trim_ltrharvest85.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_trim_ltrdigest85.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_trim_ltrdigest85");

    my @ltrh_opts = qw(-longoutput -seqids -mintsd -maxtsd -minlenltr -maxlenltr -mindistltr 
                       -maxdistltr -similar -vic -index -out -outinner -gff3);

    my @ltrh_args = ("no","yes","4","6","70","500","280","1500","85","10",
		     $index,$ltrh_out,$ltrh_out_inner,$ltrh_gff);

    @ltrh_cmd{@ltrh_opts} = @ltrh_args;
    
    my $ltr_succ  = $self->run_ltrharvest(\%ltrh_cmd);
    my $gffh_sort = $self->sort_gff($ltrh_gff) if -s $ltrh_gff;

    if (defined $gffh_sort && -s $gffh_sort) {
	my @ltrd_opts = qw(-trnas -hmms -aliout -aaout -seqfile -matchdescstart -seqnamelen -o -outfileprefix);
	my @ltrd_args = ($trnadb,$hmmdb,"yes","yes",$genome,"yes","50",$ltrg_gff,$ltrg_out);
	
	@ltrd_cmd{@ltrd_opts} = @ltrd_args;
	
	my $ltr_dig = $self->run_ltrdigest(\%ltrd_cmd, $gffh_sort);
	$self->clean_index if $self->clean;
	unlink $ltrh_gff;
	unlink $gffh_sort;
	
	return $ltrg_gff;
    }
    else {
	$self->clean_index if $self->clean;
	return 0;
    }
}


__PACKAGE__->meta->make_immutable;

1;