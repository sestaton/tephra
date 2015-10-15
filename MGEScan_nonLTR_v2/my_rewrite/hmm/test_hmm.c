#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include "hmm.h"


int main (int argc, char **argv){

  int c;
  char *hmm_file, *seq_file, *sig1_file, *sig2_file, *out_file, *phmm_dir, *out_dir;
  int debug=1;
  FILE *fp;
  HMM hmm;

  int i,j;
  int num_seqs;
  char **obs_seqs;

  char *O;
  int T;
  double **alpha, **beta;
  double pprob;
  int *vpath;

  /* read command line argument */
  while ((c=getopt(argc, argv, "m:s:r:a:o:p:d:")) != -1){

    switch (c){
    case 'm':
      hmm_file = optarg;   /* HMM */
      break;
      
    case 's':
      seq_file = optarg;   /* DNA SEQ */
      break;      
 
    case 'a':
      sig1_file = optarg;   /* APE */
      break;      
   
    case 'r':
      sig2_file = optarg;   /* RT */
      break;      

    case 'o':
      out_file = optarg;   /* output */
      break;      

    case 'p':
      phmm_dir = optarg;   /* program main dir */
      break;      

    case 'd':  
      out_dir = optarg;   /* out1 dir */ 
      break; 

    }
  }
 
  /* read initial model */
  fp = fopen(hmm_file, "r");
  get_hmm_from_file(fp, &hmm);
  fclose(fp);
  /*print_hmm(&hmm);*/
   
  /* read observation seq */
  fp = fopen(seq_file, "r");
  read_obs_seq_from_file(fp, &num_seqs, &obs_seqs);
  fclose(fp);

  /* test */
  for (i=0; i<num_seqs; i++){
    T = strlen(obs_seqs[i]);
    viterbi(&hmm, T, obs_seqs[i], &pprob, vpath, sig1_file, sig2_file, out_file, phmm_dir, out_dir); 
  }
  /* free memory */
  free_hmm(&hmm);
  
}
