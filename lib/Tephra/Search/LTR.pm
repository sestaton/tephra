package Tephra::Search::LTR;

use 5.010;
use Moo;
use Cwd;
use File::Spec;
use File::Find;
use IPC::System::Simple qw(system EXIT_ANY);
use Types::Standard     qw(Bool);
use Path::Class::File;
use Log::Any            qw($log);
use Try::Tiny;

has genome => (
      is       => 'ro',
      #isa      => 'Path::Class::File',
      required => 1,
      #coerce   => 1,
);

has hmmdb => (
      is       => 'ro',
      #isa      => 'Path::Class::File',
      required => 1,
      #coerce   => 1,
);

has trnadb => (
      is       => 'ro',
      #isa      => 'Path::Class::File',
      required => 1,
      #coerce   => 1,
);

has clean => (
      is       => 'ro',
      #isa      => Bool,
      required => 0,
      default  => sub { 1 },
);

sub ltr_search_strict {
    my $self = shift;
    my $genome = $self->genome;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $ltrh_out = $name."_ltrharvest99_pred-all";
    my $ltrh_out_inner = $name."_ltrharvest99_pred-inner";
    my $ltrh_gff = $name."_ltrharvest99.gff3";
    my $ltrg_gff = $name."_ltrdigest99.gff3";
    my $ltrg_out = $name."_ltrdigest99";

    my $index = $self->_create_ltr_index($gt, $genome);

    my $ltrh_args = "-longoutput yes -seqids yes -mintsd 4 -maxtsd 6 ";
    $ltrh_args   .= "-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 ";
    $ltrh_args   .= "-motif tgca -similar 99 -vic 10 -index $index ";
    $ltrh_args   .= "-out $ltrh_out -outinner $ltrh_out_inner -gff3 $ltrh_gff";

    my $ltr_succ  = $self->_run_ltrharvest($gt, $ltrh_args);
    my $gffh_sort = $self->_sort_gff($gt, $ltrh_gff);

    my $ltrd_args = "-trnas $trnadb -hmms $hmmdb -aliout no -aaout no -seqfile $genome ";
    $ltrd_args   .= "-matchdescstart yes -seqnamelen 50 -o $ltrg_gff ";
    $ltrd_args   .= "-outfileprefix $ltrg_out $gffh_sort";
 
    my $ltr_dig   = $self->_run_ltrdigest($gt, $ltrd_args);
}

sub ltr_search_relaxed {
    my $self = shift;
    my $genome = $self->genome;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    #my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $ltrh_out = $name."_ltrharvest85_pred-all";
    my $ltrh_out_inner = $name."_ltrharvest85_pred-inner";
    my $ltrh_gff = $name."_ltrharvest85.gff3";
    my $ltrg_gff = $name."_ltrdigest85.gff3";
    my $ltrg_out = $name."_ltrdigest85";

    my $index = $self->_create_ltr_index($gt, $genome);

    my $ltrh_args = "-longoutput yes -seqids yes -mintsd 4 -maxtsd 6 ";
    $ltrh_args   .= "-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 ";
    $ltrh_args   .= "-similar 85 -vic 10 -index $index ";
    $ltrh_args   .= "-out $ltrh_out -outinner $ltrh_out_inner -gff3 $ltrh_gff";

    my $ltr_succ  = $self->_run_ltrharvest($gt, $ltrh_args);
    my $gffh_sort = $self->_sort_gff($gt, $ltrh_gff);

    my $ltrd_args = "-trnas $trnadb -hmms $hmmdb -aliout no -aaout no -seqfile $genome ";
    $ltrd_args   .= "-matchdescstart yes -seqnamelen 50 -o $ltrg_gff ";
    $ltrd_args   .= "-outfileprefix $ltrg_out $gffh_sort";
 
    my $ltr_dig   = $self->_run_ltrdigest($gt, $ltrd_args);
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

    try {
        system([0..5], "$gt ltrharvest", $args);
    }
    catch {
        $log->error("LTRharvest failed. Here is the exception: $_\nExiting.");
        exit(1);
    };

}

sub _run_ltrdigest {
    my $self = shift;
    my ($gt, $args) = @_;

    try {
        system([0..5], "$gt ltrdigest", $args);
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
    my $gff_sort = $name."_".$suffix;
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
    find( sub { push @files, $File::Find::name if /\.tis|\.suf|\.lcp|\.ssp|\.sds|\.des|\.dna/ }, $dir);
    unlink @files;
}

1;
