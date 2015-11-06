#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>
#include <string.h>
#include <values.h>
#include <errno.h>
#include "hmm.h"
#include "util_lib.h"

int ape_state[NO_APE] = {APE_JOC, APE_I, APE_CR1, APE_Rex, APE_R1, APE_Tad1, APE_RTE, APE_L1, APE_Rand1, APE_L2};    
int id_state[NO_APE] = {ID_JOC, ID_I, ID_CR1, ID_Rex, ID_R1, ID_Tad1, ID_RTE, ID_L1, ID_Rand1, ID_L2};
int rt_state[NO_RT] = {RT_JOC, RT_I, RT_CR1, RT_Rex, RT_R1, RT_Tad1, RT_RTE, RT_L1, RT_Rand1, RT_L2, RT_CRE, RT_R2};
char *phmm_dir;
char *out_dir;

void get_hmm_from_file(FILE *fp, HMM *hmm_ptr){

  int i, j, num_trans, num_state, from, to;
  double pr;

  fscanf(fp, "Symbol= %d\n", &(hmm_ptr->M));
  fscanf(fp, "State= %d\n", &(hmm_ptr->N));

  /* default transition */
  hmm_ptr->A = (double **)dmatrix(hmm_ptr->N, hmm_ptr->N);
  for (i=0; i<hmm_ptr->N; i++){
    for (j=0; j<hmm_ptr->N; j++){
      hmm_ptr->A[i][j] = 0;
    }
  }

  /* transition */
  fscanf(fp, "Transition= %d\n", &num_trans);
  for (i=0; i<num_trans; i++){
    fscanf(fp, "%d %d %lf\n", &from, &to, &pr);
    hmm_ptr->A[from][to] = pr;
  }

  /* start state*/
  fscanf(fp, "Pi= %d\n", &num_state);
  hmm_ptr->pi = (double *)dvector(hmm_ptr->N);
  for (i=0; i<hmm_ptr->N; i++){
    fscanf(fp, "%lf\n", &pr);
    hmm_ptr->pi[i] = pr;
  }
}
 

void read_obs_seq_from_file(FILE *fp, int *num_seqs_ptr, char ***obs_seqs_ptr){

  int i,j;
  char mystring[100];

  *num_seqs_ptr = 1;          /* for now, one seq in a file */
  *obs_seqs_ptr = (char **)malloc(*num_seqs_ptr * sizeof(int));

  char **obs_seqs = *obs_seqs_ptr;
  i=0;
  j=0;
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
    }else{
      i += (strlen(mystring)-1);
    }
  }
  rewind(fp);
  obs_seqs[0] = (char *)malloc(i * sizeof(char)+1);
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
    }else{
      memcpy(obs_seqs[0]+j, mystring, strlen(mystring)-1);
      j += (strlen(mystring)-1);
    }
  }
}
  

void viterbi(HMM *hmm_ptr, int T, char *O, double *pprob, int *vpath, char *signal_file1, char *signal_file2, char *out_file, char *phmm_dir1, char *out_dir1){

  /*
    chr *O: sequecne
    int *vpath: optimal path after backtracking
    int T: the length of sequence


    hmm.A[][]: transition prob
    hmm.B[][]: emission prob
    hmm.N: the number of state
    
    int i: index for start state, current state
    int j: index for end start
    int t: index for string
    int t1: subindex for string

    int state_end;
    double state_score;
  */

  int i, j, t, t1, d, ss, si;
  double max_val, temp_val;
  int max_state;
  int **path;           /* viterbi path */
  double **alpha;       /* viterbi array */
  int *signal;          /* chromosomal position for the index in X axis of array */
  int *path_signal;     /* 0: start of fragment, 1: end */
 
  int temp1, temp2;
  double temp3, temp4;
  FILE *fp;
  int state_end;
  double state_score;
  int *ape, *rt;
  int *start, *end;
  int flag;
  char state_flag[10];
  int prev_state;
  int final_t;

  int sig_id = 0;      /* index for signal of fragments */
  int num_sig = 0;     /* total numer of signals */
  int temp_id;
  double temp_prob;
  int check_prob;
  int debug = 0;

  if ((phmm_dir = (char *)malloc(50 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of phmm_dir\n");
    exit;
  }
  if ((out_dir = (char *)malloc(50 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of out_dir\n");
    exit;
  }

  memset(state_flag,0,10);

  memset(phmm_dir,0,50);
  phmm_dir = phmm_dir1;

  memset(out_dir,0,50);
  out_dir = out_dir1; 

  /***************************************************************/
  /* initialize viterbi array                                    */
  /***************************************************************/
  ape = (int *)ivector(T+1);
  rt = (int *)ivector(T+1);
  start = (int *)ivector(hmm_ptr->N);
  end = (int *)ivector(hmm_ptr->N);
  /*****************************************************************/
  /* mark signal                                                   */
  /*****************************************************************/  
  for (t=0; t<T; t++){
    ape[t] = 0;
    rt[t] = 0;
  }
  fp = fopen(signal_file1, "r");
  while (fscanf(fp, "%d %d %lf %lf\n", &temp1, &temp2, &temp3, &temp4) != EOF){
    ape[temp1] = 1;
    ape[temp2] = 2;
    num_sig += 2;
  }
  fclose(fp);
  
  fp = fopen(signal_file2, "r");
  while (fscanf(fp, "%d %d %lf %lf\n", &temp1, &temp2, &temp3, &temp4) != EOF){
    rt[temp1] = 1;
    rt[temp2] = 2;
    num_sig += 4; 
  }
  fclose(fp);

  num_sig *= 2; 

  /***************************************************************/
  /* initialize viterbi array                                    */
  /***************************************************************/
  alpha = (double **)dmatrix(hmm_ptr->N, num_sig);
  path = (int **)imatrix(hmm_ptr->N, num_sig);
  signal = (int *)ivector(num_sig);
  path_signal = (int *)ivector(num_sig);

  /******************************************************************/
  /* initialize the first IE state                                  */
  /******************************************************************/
  start[IE]=0;
  signal[sig_id] = 0;
  sig_id ++;

  

  for (i=0; i<hmm_ptr->N; i++){
    alpha[i][0] = -1 * log(hmm_ptr->pi[i]);
    for (t=1; t<num_sig; t++){
      alpha[i][t] = DBL_MAX;
    }
  }

  /******************************************************************/
  /* do recursion to fill out rest of the columns                   */
  /******************************************************************/

  for (t = 1; t < T - 2; t++) {

    flag=0;

    /*********************/
    /* state APE start  */
    /*********************/      
    if (ape[t]==1){

      flag=1;      
      strcpy(state_flag, "APE start");
      path_signal[sig_id] = 0;

      alpha[IE][sig_id] = DBL_MAX;
      path[IE][sig_id] = IE;
      start[IE] = sig_id;

      for(ss=0; ss<NO_APE; ss++){
	alpha[ape_state[ss]][sig_id] = alpha[IE][sig_id-1] - log(hmm_ptr->A[IE][ape_state[ss]]);
	path[ape_state[ss]][sig_id] = IE;
	start[ape_state[ss]] = sig_id;
      }	
    }
    /*********************/
    /* state APE end     */
    /*********************/      
          
    if (ape[t] == 2) {

      flag=1;
      strcpy(state_flag,  "APE end");
      path_signal[sig_id] = 1;

      alpha[IE][sig_id] = DBL_MAX;
      path[IE][sig_id] = IE;

      for (ss=0; ss<NO_APE; ss++){
	state_score = get_prob(ape_state[ss], signal[start[ape_state[ss]]], t, O);
	if (state_score == DBL_MAX){
	  alpha[ape_state[ss]][sig_id] = DBL_MAX;
	}else{
	  alpha[ape_state[ss]][sig_id] = alpha[ape_state[ss]][start[ape_state[ss]]] - state_score - log(dist_ape(t-signal[start[ape_state[ss]]]+1)); 
	}
	path[ape_state[ss]][sig_id] = ape_state[ss];
      }
    }
    /*********************/
    /* state ID/IE start */
    /*********************/      
    
    if (ape[t-1] == 2){

      flag=1;
      strcpy(state_flag,  "ID/IE start");
      path_signal[sig_id] = 0;

      alpha[IE][sig_id] = alpha[IE][sig_id-1];
      path[IE][sig_id] = IE;
      start[IE] = sig_id;

      for (si=0; si<NO_APE; si++){
	temp_val = alpha[ape_state[si]][sig_id-1] - log(hmm_ptr->A[ape_state[si]][IE]);
	if (temp_val < alpha[IE][sig_id]){
	  alpha[IE][sig_id] = temp_val;
	  path[IE][sig_id] = ape_state[si];
	}
      }

      for (si=0; si<NO_APE; si++){
	if (alpha[ape_state[si]][sig_id-1]==DBL_MAX){
	  alpha[id_state[si]][sig_id]  = DBL_MAX;
	}else{
	  alpha[id_state[si]][sig_id] = alpha[ape_state[si]][sig_id-1] - log(hmm_ptr->A[ape_state[si]][id_state[si]]);
	}
	path[id_state[si]][sig_id] = ape_state[si];
	start[id_state[si]] = sig_id;
      }
    }
    /*********************/
    /* state ID/IE end   */
    /*********************/      
        
    if (rt[t+1] == 1 ) {

      flag=1;
      strcpy(state_flag,  "ID/IE end");
      path_signal[sig_id] = 1;
 
      state_score = get_prob(IE, signal[start[IE]], t, O);

      alpha[IE][sig_id] = alpha[IE][start[IE]] - state_score - log(dist_ie(t-signal[start[IE]]+1));   
      path[IE][sig_id] = IE;

      if (start[id_state[0]] >-1){
	for (si=0; si<NO_APE; si++){
	  state_score = get_prob(id_state[si], signal[start[id_state[si]]], t, O);
	  
	  if (dist_id(t-signal[start[id_state[si]]]+1)==0){
	    alpha[id_state[si]][sig_id] = DBL_MAX;   
	  }else{
	    alpha[id_state[si]][sig_id] = alpha[id_state[si]][start[id_state[si]]] - state_score - log(dist_id(t-signal[start[id_state[si]]]+1));   
	  }
	  path[id_state[si]][sig_id] = id_state[si];
	  start[id_state[si]]=-1;
	}
      }else{
	for (si=0; si<NO_APE; si++){
	  alpha[id_state[si]][sig_id] = DBL_MAX;   
	  path[id_state[si]][sig_id] = -1;
	}
      }
    }
    /*********************/
    /* state RT1 start   */
    /*********************/
           
    if (rt[t] == 1) {         

      flag=1;
      strcpy(state_flag,  "RT start");
      path_signal[sig_id] = 0;

      alpha[IE][sig_id] = DBL_MAX;
      path[IE][sig_id] = IE;
      start[IE] = sig_id;
      
      for (si=0; si<NO_RT; si++){

	alpha[rt_state[si]][sig_id] = alpha[IE][sig_id-1] - log(hmm_ptr->A[IE][rt_state[si]]);
	path[rt_state[si]][sig_id] = IE;
	start[rt_state[si]] = sig_id;	

	temp_val = alpha[id_state[si]][sig_id-1] - log(hmm_ptr->A[id_state[si]][rt_state[si]]);

	if (temp_val < alpha[rt_state[si]][sig_id]){
	  alpha[rt_state[si]][sig_id] = temp_val;
	  path[rt_state[si]][sig_id] = id_state[si];
	}
      }
    }
    /*********************/
    /* state RT1 end     */
    /*********************/      
        
    if (rt[t] == 2) {

      flag=1;
      strcpy(state_flag,  "RT end");
      path_signal[sig_id] = 1;

      alpha[IE][sig_id] = DBL_MAX;
      path[IE][sig_id] = IE;
      
      for (si=0; si<NO_RT; si++) {
	state_score = get_prob(rt_state[si], signal[start[rt_state[si]]], t, O);
	if (state_score == DBL_MAX){
	  alpha[rt_state[si]][sig_id] = DBL_MAX;
	}else{
	  alpha[rt_state[si]][sig_id] = alpha[rt_state[si]][start[rt_state[si]]] - state_score - log(dist_rt(t-signal[start[rt_state[si]]]+1));   
	}
	path[rt_state[si]][sig_id] = rt_state[si];
      }
    }
    
    /*********************/
    /* state IE start    */
    /*********************/      
        
    if (rt[t-1] == 2){

      flag=1;
      strcpy(state_flag,  "IE start");
      path_signal[sig_id] = 0;

      alpha[IE][sig_id] = DBL_MAX;
      path[IE][sig_id] = IE;
      start[IE] = sig_id;

      for (si=0; si<NO_RT; si++){
	temp_val = alpha[rt_state[si]][sig_id-1] - log(hmm_ptr->A[rt_state[si]][IE]);
	if (temp_val < alpha[IE][sig_id]){
	  alpha[IE][sig_id] = temp_val;
	  path[IE][sig_id] = rt_state[si];
	}
      }
      for (si=0; si<NO_APE; si++){
	alpha[id_state[si]][sig_id]=DBL_MAX;
	path[id_state[si]][sig_id] = -1;
      }
    }
    /*********************/
    /* state IE end     */
    /*********************/      
        
    if (ape[t+1] == 1  ) {

      flag=1;
      strcpy(state_flag,  "IE end");
      path_signal[sig_id] = 1;

      state_score = get_prob(IE, signal[start[IE]], t, O);

      if (dist_ie(t-signal[start[IE]]+1)==0){
	alpha[IE][sig_id] = DBL_MAX;
      }else{
	alpha[IE][sig_id] = alpha[IE][start[IE]] - state_score - log(dist_ie(t-signal[start[IE]]+1));
      }
      path[IE][sig_id] = IE;
    }

    /************************/
    /* print                */
    /************************/
    
    if (flag==1){
      final_t = sig_id;
      signal[sig_id] = t;
      if (debug==1){
	printf("%d %d %s\t\t", sig_id, t, state_flag);
	
	for (i=0; i<hmm_ptr->N; i++){
	  if (alpha[i][sig_id] < DBL_MAX ){
	    printf("%lf\t", alpha[i][sig_id]);
	  }else{
	    printf("-\t");
	  }
	}
	printf("\n\n");

	printf("%d %d %s\t\t", sig_id, t, state_flag);
	
	for (i=0; i<hmm_ptr->N; i++){
	  printf("%d\t", path[i][sig_id]);
	}
	printf("\n\n\n");
      }
      sig_id++;
    }

  }

  /***********************************************************/
  /* backtrack array to find the optimal path                */
  /***********************************************************/
  
  *pprob = DBL_MAX;
  vpath = (int *)ivector(num_sig);

  for (i = 0; i < hmm_ptr->N; i++){
    if (alpha[i][final_t] < *pprob && alpha[i][final_t] > 0.00001){
      *pprob = alpha[i][final_t];
      vpath[final_t] = i;
    }
  } 
  if (*pprob == DBL_MAX) {
    vpath[final_t]=0;
  }

  fp = fopen (out_file , "w");
  if (fp == NULL) {
    fprintf(stderr, "Can't open input file out_file\n");
    exit(1);
  }
  for(sig_id=final_t-1; sig_id>=0; sig_id--){    
    vpath[sig_id] = path[vpath[sig_id+1]][sig_id+1];
 
    if (path_signal[sig_id] == 1){
      fprintf(fp, "%d %d %lf\n", signal[sig_id], vpath[sig_id], alpha[vpath[sig_id]][sig_id]);
    }
  }
  fclose(fp);
  //if (ret >
 
  free_ivector(ape);
  free_ivector(rt);
  if (start) {
    free_ivector(start);
  }
  if (end) {
    free_ivector(end);
  }
  if (hmm_ptr->N) {
    if (path) {
      free_imatrix(path, hmm_ptr->N);
    }
    if (alpha) {
      free_dmatrix(alpha, hmm_ptr->N);
    }
  }
  if (signal) {
    free_ivector(signal);
  }
  if (path_signal) {
    free_ivector(path_signal);
  }
}

double dist_ape(int dist){

  double result;

  if (400 <= dist && dist < 500){
    result = 0.01;
  }else if (500 <= dist && dist < 600){
    result = 0.01;
  }else if (600 <= dist && dist < 700){
    result = 0.77;
  }else if (700 <= dist && dist < 800){
    result = 0.20;
  }else{
    result = 0.0000000000000001;
  }
  return result;
}


double dist_id(int dist){

  double result;

  if (300<=dist && dist<600){
    result = 0.12;
  }else if (600<= dist && dist< 700){
    result = 0.84;
  }else if (700<=dist && dist<800){
    result = 0.02;
  }else if (800<= dist && dist< 900){
    result = 0.01;
  }else if (900<= dist && dist< 1000){
    result = 0.01;
  }else if (1000 <= dist && dist< 1100){
    result = 0.01;
  }else{
    result = 0.000000000000000001;
  }
  return result;
}

double dist_rt(int dist){

  double result;

  if (700 <= dist && dist < 800){
    result = 0.01;
  }else if (800 <= dist && dist < 900){
    result = 0.02;
  }else if (900 <= dist && dist < 1000){
    result = 0.01;
  }else if (1000 <= dist && dist < 1100){
    result = 0.01;
  }else if (1100 <= dist && dist < 1200){
    result = 0.04;
  }else if (1200 <= dist && dist < 1300){
    result = 0.74;
  }else if (1300 <= dist && dist < 1400){
    result = 0.17;
  }else if (1400 <= dist && dist < 1500){
    result = 0.01;
  }else{
    result = 0.0000000000000000001;
  }
  return result;
}

double dist_ie(int dist){

  double result;

  if (1500 <= dist && dist <3000){
    result = 0.15;
  }else if (1 <= dist && dist < 100){
    result = 0.14;
  }else if (100 <= dist && dist < 1500){
    result = 0.01;
  }else{
    result=0.7;
  }
  return result;
}


double  get_prob(int state, int start_pos, int end_pos,  char *O){

  double state_score=0;
  int seq_len;
  int flag = 0;   /* 1: ape or rt  2: id */
  int i;

  for (i=0; i<NO_APE; i++){
    if (state == ape_state[i]){
      flag=1;
    }
  }
  for (i=0; i<NO_RT; i++){
    if (state == rt_state[i]){
      flag=1;
    }
  }
  for (i=0; i<NO_APE; i++){
    if (state == id_state[i]){
      flag=2;
    }
  }

  if (state == IE){

    // long fragments should be IE not ID
    if (end_pos-start_pos > 10000){
      state_score = 0.00000000001;   
    }else{
      get_hydro(start_pos, end_pos, &state_score, O);
    }
    state_score = log(1-state_score); 

  }else if (flag == 1){

    get_phmm(state, start_pos, end_pos, &state_score, O, &seq_len);

    if (state_score == -1){
      state_score = DBL_MAX;
    }else{
      state_score = log(state_score);
    }

  }else if (flag == 2){

    if (end_pos-start_pos > 10000){
      state_score = 0.00000000001;
    }else{
      get_hydro(start_pos, end_pos, &state_score, O);
    }
    state_score = log(state_score); 
  }

  return state_score;
}

void get_hydro(int start, int end, double *score, char *O){


  double hydro_kd[20] = {4.5, 3.9, 3.5, 3.5, 3.5, 3.5, 3.2, 1.6, 1.3, 0.9, 
			 0.8, 0.7, 0.4, -1.8, -1.9, -2.5, -2.8, -3.8, -4.2, -4.5};
  double hydro_ww[20] = {0.81, 0.99, 0.42, 1.23, 0.58, 2.02, 0.96, 0.45, -0.94, -1.85, 
			 0.13, 0.14, 0.01, 0.17, -0.23, -0.24, -1.13, -0.56, 0.07, -0.31};
  double hydro_hh[20] = {2.58, 2.71, 2.05, 3.49, 2.36, 2.68, 2.06, 2.23, 0.68, 0.3, 
			 0.84, 0.52, 0.74, 0.11, -0.1, -0.13, -0.32, -0.55, -0.31, -0.6};

  double score_kd, score_ww, score_hh, temp_kd, temp_ww, temp_hh;
  double h_kd, h_ww, h_hh;
  double r_kd, r_ww, r_hh;
  double p_kd, p_ww, p_hh;

  int count_stop=10000, count_aa, temp_stop, temp_aa;
  char temp_file1[100];
  char temp_file2[100];
  char *command;
  char *domain_bp;
  char *st;
  char *seq;
  int i;
  FILE *fp;
  char buffer[10]; /* long enough to hold your number + null */

  if ((command = (char *)malloc(10000 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of command\n");
    exit;
  }
  memset(command, '\0', 10000);

  if ((domain_bp = (char *)malloc(10000 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of domain_bp\n");
    exit;
  }
  memset(domain_bp, '\0', 10000);

  if ((st = (char *)malloc(70 * sizeof(char)))==NULL){
    printf("ERROR: allocation of st\n");
    exit;
  }
  memset(st, '\0', 70);

  if ((seq = (char *)malloc(20000 * sizeof(char)))==NULL){
    printf("ERROR: allocation of seq\n");
    exit;
  }
  memset(seq, '\0', 20000);

  strcpy(temp_file1, out_dir);
  strcat(temp_file1, "ppppp");

  strcpy(temp_file2, out_dir);
  strcat(temp_file2, "qqqqq");

  if (end-start+5 > 10000){
    end = start + 9995;
  }

  memcpy(domain_bp, &O[start], end-start+1);

  if ((fp = fopen(temp_file1, "w"))==NULL){
    printf("ERROR: file opening for temp_file1: %s\n", strerror(errno));
  }
  fprintf(fp, "%s\n", domain_bp);
  if(ferror(fp)) {
    printf("ERROR: file writing for domain_bp\n"); 
  }
  fclose(fp);
  if(ferror(fp)) {
    printf("ERROR: file writing for domain_bp\n");
  }

  strcpy(command, "transeq -frame=f ");
  strcat(command, temp_file1);
  strcat(command, " -outseq=");
  strcat(command, temp_file2);
  strcat(command, " 2>/dev/null");
  system(command);
  

  temp_aa = 0;
  temp_stop = 0;
  temp_kd = 0;
  temp_ww = 0;
  temp_hh = 0;

  fp = fopen(temp_file2, "r");
  while (fscanf(fp, "%s\n", st) != EOF){

    if (st[2] == '1'){

    }else if (st[2] == '2' || st[2] == '3'){

      /* calculate hydrophobic index for frame1 and frame2 */
      for ( i=0; i<strlen(seq); i++){

	if (seq[i] == 'I' || seq[i] == 'V' || seq[i] == 'L' || seq[i] == 'F' ||
	    seq[i] == 'C' || seq[i] == 'M' || seq[i] == 'A' || seq[i] == 'G' ||
	    seq[i] == 'T' || seq[i] == 'S' || seq[i] == 'W' || seq[i] == 'Y' ||
	    seq[i] == 'P' || seq[i] == 'H' || seq[i] == 'E' || seq[i] == 'Q' ||
	    seq[i] == 'D' || seq[i] == 'N' || seq[i] == 'K' || seq[i] == 'R' ) {
	  temp_kd += hydro_kd[convert_aa_int(seq[i])];
	  temp_ww += hydro_ww[convert_aa_int(seq[i])];
	  temp_hh += hydro_hh[convert_aa_int(seq[i])];
	  temp_aa ++;

	}else if  (seq[i] == '*'){
	  temp_stop++;
	}
      }

      if (temp_stop < count_stop ){
	score_kd = temp_kd;
	score_ww = temp_ww;
	score_hh = temp_hh;
	count_aa = temp_aa;
	count_stop = temp_stop;
      }
      
      temp_aa = 0;
      temp_stop = 0;
      temp_kd = 0;
      temp_ww = 0;
      temp_hh = 0;
      memset(seq, '\0', 10000);

    }else{
      strcat(seq, st);
    }
  }
  fclose(fp);

  /* the last one */
  for (i=0; i<strlen(seq); i++){
    
    if (seq[i] == 'I' || seq[i] == 'V' || seq[i] == 'L' || seq[i] == 'F' ||
	seq[i] == 'C' || seq[i] == 'M' || seq[i] == 'A' || seq[i] == 'G' ||
	seq[i] == 'T' || seq[i] == 'S' || seq[i] == 'W' || seq[i] == 'Y' ||
	seq[i] == 'P' || seq[i] == 'H' || seq[i] == 'E' || seq[i] == 'Q' ||
	seq[i] == 'D' || seq[i] == 'N' || seq[i] == 'K' || seq[i] == 'R' ) {
      temp_kd += hydro_kd[convert_aa_int(seq[i])];
      temp_ww += hydro_ww[convert_aa_int(seq[i])];
      temp_hh += hydro_hh[convert_aa_int(seq[i])];
      temp_aa ++;
      
    }else if  (seq[i] == '*'){
      temp_stop++;
    }
  }
  
  if (temp_stop < count_stop ){
    score_kd = temp_kd;
    score_ww = temp_ww;
    score_hh = temp_hh;
    count_aa = temp_aa;
    count_stop = temp_stop;
  }
  
  /* calculate the prob */
  if (count_aa>0){
    score_kd /= count_aa;
    score_ww /= count_aa;
    score_hh /= count_aa;
  }else{
    score_kd = -4.5;
    score_ww = -1.85;
    score_hh = -0.6;
  }

  h_kd = A_H_KD * exp(-1*pow(score_kd-MU_H_KD,2)/(2*pow(SIGMA_H_KD,2)));
  h_ww = A_H_WW * exp(-1*pow(score_ww-MU_H_WW,2)/(2*pow(SIGMA_H_WW,2)));
  h_hh = A_H_HH * exp(-1*pow(score_hh-MU_H_HH,2)/(2*pow(SIGMA_H_HH,2)));
  r_kd = A_R_KD * exp(-1*pow(score_kd-MU_R_KD,2)/(2*pow(SIGMA_R_KD,2)));
  r_ww = A_R_WW * exp(-1*pow(score_ww-MU_R_WW,2)/(2*pow(SIGMA_R_WW,2)));
  r_hh = A_R_HH * exp(-1*pow(score_hh-MU_R_HH,2)/(2*pow(SIGMA_R_HH,2)));
  
  p_kd = h_kd / (h_kd + r_kd);
  p_ww = h_ww / (h_ww + r_ww);
  p_hh = h_hh / (h_hh + r_hh);

  if (p_kd>p_ww){
    *score = p_kd;
  }else{
    *score = p_ww;
  }
  if (p_hh>*score){
    *score = p_hh;
  }

  if (count_aa==0){
    *score = 0.000000000000001;
  }

   free(st);
   free(seq);
   free(command);
   free(domain_bp);
}

void get_phmm(int state, int start, int end, double *score, char *O, int *seq_len){

  char *domain_bp;
  char *command;
  char *hmm_file;
  FILE *fp;
  char temp[10];
  //char *temp = malloc(10 * sizeof(char));

  if ((command = (char *)malloc(10000 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of command\n");
    exit;
  }
  if ((domain_bp = (char *)malloc(10000 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of domain_bp\n");
    exit;
  }
  if ((hmm_file = (char *)malloc(100 * sizeof(char)))==NULL){ 
    printf("ERROR: allocation of hmm_file\n");
    exit;
  }
  memset(command, '\0', 10000);
  if (end-start+5 > 10000){
    end = start + 9995;
  }
  memset(domain_bp, '\0', 10000);
  memcpy(domain_bp, &O[start], end-start+1);

  strcpy(hmm_file, phmm_dir);
  strcat(hmm_file, "pHMM/");
  if (state == RT_JOC){    
    strcat(hmm_file, "Jockey.rt.hmm ");
  }else if (state == APE_JOC){
    strcat(hmm_file, "Jockey.en.hmm ");
  }else if (state == RT_I){
    strcat(hmm_file, "I.rt.hmm ");
  }else if (state == APE_I){
    strcat(hmm_file, "I.en.hmm ");
  }else if (state == RT_CR1){
    strcat(hmm_file, "CR1.rt.hmm ");
  }else if (state == APE_CR1){
    strcat(hmm_file, "CR1.en.hmm ");
  }else if (state == RT_R1){
    strcat(hmm_file, "R1.rt.hmm ");
  }else if (state == APE_R1){
    strcat(hmm_file, "R1.en.hmm ");
  }else if (state == RT_Tad1){
    strcat(hmm_file, "Tad1.rt.hmm ");
  }else if (state == APE_Tad1){
    strcat(hmm_file, "Tad1.en.hmm ");
  }else if (state == RT_RTE){
    strcat(hmm_file, "RTE.rt.hmm ");
  }else if (state == APE_RTE){
    strcat(hmm_file, "RTE.en.hmm ");
  }else if (state == RT_L1){
    strcat(hmm_file, "L1.rt.hmm ");
  }else if (state == APE_L1){
    strcat(hmm_file, "L1.en.hmm ");
  }else if (state == RT_Rand1){
    strcat(hmm_file, "RandI.rt.hmm ");
  }else if (state == APE_Rand1){
    strcat(hmm_file, "RandI.en.hmm ");
  }else if (state == RT_L2){
    strcat(hmm_file, "L2.rt.hmm ");
  }else if (state == APE_L2){
    strcat(hmm_file, "L2.en.hmm ");
  }else if (state == RT_Rex){
    strcat(hmm_file, "Rex.rt.hmm ");
  }else if (state == APE_Rex){
    strcat(hmm_file, "Rex.en.hmm ");
  }else if (state == RT_CRE){
    strcat(hmm_file, "CRE.rt.hmm ");
  }else if (state == RT_R2){
    strcat(hmm_file, "R2.rt.hmm ");
  }

  //strcpy(command, phmm_dir);
  strcat(command, "tephra-getphmm -s ");
  strcat(command, domain_bp);
  strcat(command, " -h ");
  strcat(command, hmm_file);
  strcat(command, " -o ");
  strcat(command, out_dir);

  fp = popen(command, "r"); 
  fscanf(fp, "%s", temp);
  *score = atof(temp);
  if (fp) {
    pclose(fp);
  }
  //free(temp);

  if (*score >=0.1){
    *score=-1;
  }else if (*score == 0.0){
    *score = log10(1.0e-299)/-300;
  }else{
        *score = log10(*score)/-300;
  }

  if (domain_bp) {
    free(domain_bp);
  }
  if (command) {
    free(command);
  }
  if (hmm_file) {
    free(hmm_file);
  }

}
 
void free_hmm(HMM *hmm_ptr){

  free_dmatrix(hmm_ptr->A, hmm_ptr->N);
  free_dmatrix(hmm_ptr->B, hmm_ptr->N);
  free_dvector(hmm_ptr->pi);
}

