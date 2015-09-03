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

with 'Tephra::Role::GT';

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has hmmdb => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 0,
      coerce   => 1,
);

has trnadb => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 0,
      coerce   => 1,
);

has clean => (
      is       => 'ro',
      isa      => 'Bool',
      required => 0,
      default  => 1,
);

sub ltr_search_strict {
    my $self = shift;
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;
    
    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_ltrharvest99_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_ltrharvest99_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_ltrharvest99.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_ltrdigest99.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_ltrdigest99");

    my $index = $self->_create_ltr_index($gt, $genome);

    my $ltrh_args = "-longoutput no -seqids yes -tabout no -mintsd 4 -maxtsd 6 ";
    $ltrh_args   .= "-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 ";
    $ltrh_args   .= "-motif tgca -similar 99 -vic 10 -index $index ";
    $ltrh_args   .= "-out $ltrh_out -outinner $ltrh_out_inner -gff3 $ltrh_gff";

    my $ltr_succ  = $self->_run_ltrharvest($gt, $ltrh_args);
    my $gffh_sort = $self->_sort_gff($gt, $ltrh_gff);

    my $ltrd_args = "-trnas $trnadb -hmms $hmmdb -aliout no -aaout no -seqfile $genome ";
    $ltrd_args   .= "-matchdescstart yes -seqnamelen 50 -o $ltrg_gff ";
    $ltrd_args   .= "-outfileprefix $ltrg_out $gffh_sort";
 
    my $ltr_dig   = $self->_run_ltrdigest($gt, $ltrd_args);
    $self->_clean_index if $self->clean;
    unlink $ltrh_gff;
    unlink $gffh_sort;
    
    return $ltrg_gff;
}

sub ltr_search_relaxed {
    my $self = shift;
    my $genome = $self->genome->absolute;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    if ($name =~ /(\.fa.*)/) {
	$name =~ s/$1//;
    }
    
    my $ltrh_out = File::Spec->catfile($path, $name."_ltrharvest85_pred-all");
    my $ltrh_out_inner = File::Spec->catfile($path, $name."_ltrharvest85_pred-inner");
    my $ltrh_gff = File::Spec->catfile($path, $name."_ltrharvest85.gff3");
    my $ltrg_gff = File::Spec->catfile($path, $name."_ltrdigest85.gff3");
    my $ltrg_out = File::Spec->catfile($path, $name."_ltrdigest85");

    my $index = $self->_create_ltr_index($gt, $genome);

    my $ltrh_args = "-longoutput no -seqids yes -tabout no -mintsd 4 -maxtsd 6 ";
    $ltrh_args   .= "-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 ";
    $ltrh_args   .= "-similar 85 -vic 10 -index $index ";
    $ltrh_args   .= "-out $ltrh_out -outinner $ltrh_out_inner -gff3 $ltrh_gff";

    my $ltr_succ  = $self->_run_ltrharvest($gt, $ltrh_args);
    my $gffh_sort = $self->_sort_gff($gt, $ltrh_gff);

    my $ltrd_args = "-trnas $trnadb -hmms $hmmdb -aliout no -aaout no -seqfile $genome ";
    $ltrd_args   .= "-matchdescstart yes -seqnamelen 50 -o $ltrg_gff ";
    $ltrd_args   .= "-outfileprefix $ltrg_out $gffh_sort";
 
    my $ltr_dig   = $self->_run_ltrdigest($gt, $ltrd_args);
    $self->_clean_index if $self->clean;
    unlink $ltrh_gff;
    unlink $gffh_sort;

    return $ltrg_gff;
}

sub _create_ltr_index {
    my $self = shift;
    my ($gt, $genome) = @_;

    #my $gt = $self->get_gt_exec;
    my $index = $genome.".index";
    my $index_cmd = "$gt suffixerator ";
    $index_cmd .= "-db $genome -indexname $index -tis -suf -lcp -ssp -sds -des -dna";
    try {
	system([0..5], $index_cmd);
    }
    catch {
	$log->error("Unable to make suffixerator index. Here is the exception: $_\nExiting.");
        exit(1);
    };

    return $index;
}

sub _run_ltrharvest {
    my $self = shift;
    my ($gt, $args) = @_;

    my $ltrh_cmd = "$gt ltrharvest $args 2>&1 > /dev/null";
    #say STDERR $ltrh_cmd;

    try {
        system([0..5], $ltrh_cmd);
    }
    catch {
        $log->error("LTRharvest failed. Here is the exception: $_\nExiting.");
        exit(1);
    };

}

sub _run_ltrdigest {
    my $self = shift;
    my ($gt, $args) = @_;

    my $ltrd_cmd = "$gt ltrdigest $args";
    try {
        system([0..5], $ltrd_cmd);
    }
    catch {
        $log->error("LTRdigest failed. Here is the exception: $_\nExiting.");
        exit(1);
    };

}

sub _sort_gff {
    my $self = shift;
    my ($gt, $gff) = @_;

    my ($name, $path, $suffix) = fileparse($gff, qr/\.[^.]*/);
    my $gff_sort = File::Spec->catfile($path, $name."_sort.".$suffix);
    #my $gt = $self->get_gt_exec;

    my $sort_cmd = "$gt gff3 -sort $gff > $gff_sort";
    try {
        system([0..5], $sort_cmd);
    }
    catch {
        $log->error("'gt gff3 -sort' failed. Here is the exception: $_\nExiting.");
        exit(1);
    };

    return $gff_sort;
}

sub _clean_index {
    my $self = shift;

    my $dir = getcwd();
    my @files;
    find( sub { push @files, $File::Find::name 
		    if /\.llv|\.md5|\.prf|\.tis|\.suf|\.lcp|\.ssp|\.sds|\.des|\.dna|\.esq|\.prj|\.ois/ 
	  }, $dir);
    unlink @files;
}

__PACKAGE__->meta->make_immutable;

1;
