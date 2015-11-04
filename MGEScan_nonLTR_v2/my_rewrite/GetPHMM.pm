package GetPHMM;

use 5.010;
use Moose;
use MooseX::Types::Path::Class;
use File::Spec;

has outdir => ( is => 'ro', isa => 'Path::Class::Dir', required => 1, coerce => 1 );
has fasta  => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );
has phmm   => ( is => 'ro', isa => 'Path::Class::File', required => 1, coerce => 1 );

sub get_phmm {
    my $self = shift;
    my $out_dir   = $self->outdir;
    my $hmmsearch = $self->_find_hmmsearch;
    my $seq       = $self->fasta;
    my $phmm_file = $self->phmm;

    #my $seq; 
    #my $phmm_file;
    #my $out_dir;
    #my $seq_file;
    #my $pep_file;
    #my $command;
    #my @hmm_results;
    #my $hmm_result;

    #GetOptions(    'seq=s' => \$seq,
    # 	       'hmmfile=s' => \$phmm_file,#
    #'odir=s' => \$out_dir#,
    #           );

    my $seq_file = File::Spec->catfile($out_dir, "aaaaa");
    my $pep_file = File::Spec->catfile($out_dir, "bbbbb");

    system("echo ".$seq." > ".$seq_file);
    system("transeq -frame=f $seq_file -outseq=$pep_file 2>/dev/null");
    my $command = "$hmmsearch $phmm_file $pep_file";
    my $hmm_result = `$command`;
    
    my @hmm_results = split /\n/, $hmm_result;
    for (my $i=0; $i<=$#hmm_results; $i++){
	
	if ($hmm_results[$i] =~ /^\-\-\-\-\-\-\-\-\s\-\-\-\-\-\-\-\s/){
	    if ($hmm_results[$i+1] =~ /^\S/){
		my @temp = split /\s+/, $hmm_results[$i+1];
		print $temp[9];
	    }
	    else {
		print "1";
	    }
	    last;
	}
    }
    
}

sub _find_hmmsearch {
    my $self = shift;

    my $hmmsearch = '/home/statonse/github/tephra/MGEScan_nonLTR_v2/hmmer-2.3.2/src/hmmsearch';
    return $hmmsearch;
}

__PACKAGE__->meta->make_immutable;

1;
