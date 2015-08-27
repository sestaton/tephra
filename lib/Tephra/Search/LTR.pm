package Tephra::Search::LTR;

use 5.010;
use Moo;
use Cwd;
use File::Spec;
use File::Find;
use IPC::System::Simple qw(system EXIT_ANY);
use Types::Standard     qw( );
use Log::Any            qw($log);
use Try::Tiny;

has genome => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has hmmdb => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has trnadb => (
      is       => 'ro',
      isa      => 'Path::Class::File',
      required => 1,
      coerce   => 1,
);

has clean => (
      is       => 'ro',
      isa      => Bool,
      required => 0,
      default  => 1,
);

sub ltr_search_strict {
    my $self = shift;
    my $genome = $self->genome;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $ltrh_gff = $name."_ltrharvest99.gff3";
    my $ltrg_gff = $name."_ltrdigest99.gff3";

    #gff_sort=${dbbase}_ltrharvest99_sort.gff3
    #gff_h=${dbbase}_ltrdigest99.gff3

    my $index = $self->_create_ltr_index($gt, $genome);

    my $ltrh_args = "-longoutput yes -seqids yes -mintsd 4 -maxtsd 6 ";
    $ltrh_args   .= "-minlenltr 100 -maxlenltr 6000 -mindistltr 1500 -maxdistltr 25000 ";
    $ltrh_args   .= "-motif tgca -similar 99 -vic 10 -index $index -out Ha412v1r1_ltrharvest99_pred-all ";
    $ltrh_args   .= "-outinner Ha412v1r1_ltrharvest99_pred-inner -gff3 $gff";

    my $ltr_succ  = $self->_run_ltrharvest($gt, $ltrh_args);
    my $gffh_sort = $self->_sort_gff($gt, $gff);

    my $ltrd_args = "-trnas $trnas -hmms $hmm -aliout yes -aaout yes -seqfile $db ";
    $ltrd_args   .= "-matchdescstart yes -seqnamelen 50 -o $gff_h ";
    $ltrd_args   .= "-outfileprefix Ha412v1r1_ltrdigest99 $gff_sort";
 
    my $ltr_dig   = $self->_run_ltrdigest($gt, $ltrd_args);

}

sub ltr_search_relaxed {
    my $self = shift;
    my $genome = $self->genome;
    my $hmmdb  = $self->hmmdb;
    my $trnadb = $self->trnadb;
    my $gt = $self->get_gt_exec;

    my ($name, $path, $suffix) = fileparse($genome, qr/\.[^.]*/);
    my $ltrh_gff = $name."_ltrharvest85.gff3";
    my $ltrg_gff = $name."_ltrdigest85.gff3";

    db=Ha412v1r1_genome_no_cp-mt-rd_chr-q.fasta
	dbbase=$(echo ${db%.*})
    gff=${dbbase}_ltrharvest85.gff3
    gff_sort=${dbbase}_ltrharvest85_sort.gff3
    gff_h=${dbbase}_ltrdigest85.gff3
    index=${db}.index
    hmm=/db/transposable+element_hmms/transposable+element-2.hmm
    trnas=/db/eukaryotic-tRNAs.fas

    #time $gt suffixerator -db $db -indexname $index -tis -suf -lcp -ssp -sds -des -dna

    # ltrharvest
    time $gt ltrharvest \
    -longoutput yes \
    -seqids yes \
    -mintsd 4 \
    -maxtsd 6 \
    -minlenltr 100 \
    -maxlenltr 6000 \
    -mindistltr 1500 \
    -maxdistltr 25000 \
    -similar 85 \
    -vic 10 \
    -index $index \
    -out Ha412v1r1_ltrharvest85_pred-all \
    -outinner Ha412v1r1_ltrharvest85_pred-inner \
    -gff3 $gff
    

    # sort the gff3 file(s) prior to running ltrdigest
    time $gt gff3 -sort $gff > $gff_sort

    # ltrdigest
    time $gt ltrdigest \
    -trnas $trnas \
    -hmms $hmm \
    -aliout yes \
    -aaout yes \
    -seqfile $db \
    -matchdescstart yes \
    -seqnamelen 50 \
    -o $gff_h \
    -outfileprefix Ha412v1r1_ltrdigest85 $gff_sort

}

sub _create_ltr_index {
    my $self = shift;
    my ($gt, $genome) = @_;

    my $gt = $self->get_gt_exec;
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
    my $gt = $self->get_gt_exec;

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
