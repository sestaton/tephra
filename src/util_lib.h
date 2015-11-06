#include <stdlib.h>
#include <stdio.h>
#include <math.h>


double **dmatrix(int nrh, int nch);
int **imatrix(int nrh, int nch);
double *dvector(int nh);
int *ivector(int nh);
void free_dvector(double *v);
void free_dmatrix(double **m,int nrh);
void free_imatrix(int **m, int nrh);
void free_ivector(int *v);
int convert_aa_int (char nt);
