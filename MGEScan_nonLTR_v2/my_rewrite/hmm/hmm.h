#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define NO_RT 12
#define NO_APE 10

#define NOSTATE -1
#define IE 0
#define APE_JOC 1
#define ID_JOC 2
#define RT_JOC 3
#define APE_I 4
#define ID_I 5
#define RT_I 6
#define APE_CR1 7
#define ID_CR1 8
#define RT_CR1 9
#define APE_Rex 10
#define ID_Rex 11
#define RT_Rex 12
#define APE_R1 13
#define ID_R1 14
#define RT_R1 15
#define APE_Tad1 16
#define ID_Tad1 17
#define RT_Tad1 18
#define APE_RTE 19
#define ID_RTE 20
#define RT_RTE 21
#define APE_L1 22
#define ID_L1 23
#define RT_L1 24
#define APE_Rand1 25
#define ID_Rand1 26
#define RT_Rand1 27
#define APE_L2 28
#define ID_L2 29
#define RT_L2 30
#define RT_CRE 31
#define RT_R2 32


#define A_H_KD 0.0703
#define A_H_WW 0.0976
#define A_H_HH 0.0639
#define MU_H_KD 0.6187
#define MU_H_WW 0.2774
#define MU_H_HH 1.1681
#define SIGMA_H_KD 0.3514
#define SIGMA_H_WW 0.0859
#define SIGMA_H_HH 0.1334
#define A_R_KD 0.0819
#define A_R_WW 0.1243
#define A_R_HH 0.0794
#define MU_R_KD 0.0852
#define MU_R_WW 0.1201
#define MU_R_HH 0.9087
#define SIGMA_R_KD 0.2350
#define SIGMA_R_WW 0.0608
#define SIGMA_R_HH 0.0969

#define CODE 1



typedef struct {
  int N;          /* number of states;  Q={1,2,...,N} */
  int M;          /* number of observation symbols; V={1,2,...,M}*/
  double  **A;    /* A[1..N][1..N]. a[i][j] is the transition prob                                                                     
                           of going from state i at time t to state j                                                                        
                           at time t+1 */
  double  **B;    /* B[1..N][1..M]. b[j][k] is the probability of                                                                      
		     of observing symbol k in state j */
  double  *pi;    /* pi[1..N] pi[i] is the initial state distribution. */
} HMM;




void print_hmm(HMM *hmm_ptr);
void get_hmm_from_file(FILE *fp, HMM *hmm_ptr);
void free_hmm(HMM *hmm);
void read_obs_seq_from_file(FILE *fp, int *len_seqs_ptr, char ***obs_seqs_ptr);
void viterbi(HMM *hmm_ptr,int T, char *O,  double *pprob, int *vpath, char *file1, char *file2, char *outfile, char *phmm, char *outdir);
void get_phmm(int state, int start, int end, double *score, char *O, int *seq_len);  
void get_hydro(int start, int end, double *score, char *O);
double dist_ape(int dist);
double dist_rt(int dist);
double dist_id(int dist);
double dist_ie(int dist);
double get_prob(int state, int start, int end,  char *O);
 
