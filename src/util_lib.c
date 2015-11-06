#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <malloc.h>

double **dmatrix(int num_row, int num_col){

  int i, j;
  double **m;

  m=(double **) malloc(num_row * sizeof(double*));
  if (!m) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in dmatrix()");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<num_row; i++) {
    m[i]=(double *) malloc(num_col * sizeof(double));
    if (!m[i]) {
      fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i, " in dmatrix()");
      exit(EXIT_FAILURE);
    }

    for(j=0; j<num_col; j++){
      m[i][j] = 0.0;
    }
  }
  return m;
}


int **imatrix(int num_row, int num_col){

  int i,j;
  int **m;

  m=(int **) malloc(num_row * sizeof(int*));
  if (!m) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure for points to rows in imatrix()");
    exit(EXIT_FAILURE);
  }

  for(i=0; i<num_row; i++) {
    m[i]=(int *) malloc(num_col * sizeof(int));
    if (!m[i]) {
      fprintf(stderr, "%s %d %s\n", "ERROR: Allocation failure for the row ", i ," in imatrix()");
      exit(EXIT_FAILURE);
    }

    for(j=0; j<num_col; j++){
      m[i][j] = 0;
    }
  }
  return m;
}

double *dvector(int nh){

  int j;
  double *v;

  v=(double *)malloc(nh * sizeof(double));

  if (!v) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure in dvector()");
    exit(EXIT_FAILURE);
  }

  for(j=0; j<nh; j++){
    v[j] = 0.0;
  }
  return v;
}

int *ivector(int nh){

  int j;
  int *v;

  v=(int *)malloc(nh * sizeof(int));

  if (!v) {
    fprintf(stderr, "%s\n", "ERROR: Allocation failure in ivector()");
    exit(EXIT_FAILURE);
  }

  for(j=0; j<nh; j++){
    v[j] = 0;
  }
  return v;
}


void free_dvector(double *v){
  
  free(v);
}

void free_ivector(int *v){
  
  free(v);
}

void free_dmatrix(double **m, int num_row){

  int i;

  for(i=num_row-1; i>=0; i--)
    free(m[i]);

  free(m);
}
void free_imatrix(int **m,int num_row){

  int i;

  for(i=num_row-1; i>=0; i--)
    free(m[i]);

  free(m);
}


int convert_aa_int (char nt){

  int result;

  switch(nt){
  case 'R':
    result=0;
    break;
  case 'K':
    result=1;
    break;
  case 'N':
    result=2;
    break;
  case 'D':
    result=3;
    break;
  case 'Q':
    result=4;
    break;
  case 'E':
    result=5;
    break;
  case 'H':
    result=6;
    break;
  case 'P':
    result=7;
    break;
  case 'Y':
    result=8;
    break;
  case 'W':
    result=9;
    break;
  case 'S':
    result=10;
    break;
  case 'T':
    result=11;
    break;
  case 'G':
    result=12;
    break;
  case 'A':
    result=13;
    break;
  case 'M':
    result=14;
    break;
  case 'C':
    result=15;
    break;
  case 'F':
    result=16;
    break;
  case 'L':
    result=17;
    break;
  case 'V':
    result=18;
    break;
  case 'I':
    result=19;
    break;
  }  

  return result;
}




