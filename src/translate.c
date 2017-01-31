#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>


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


void mysubstr(char *dest, char *src, int position, int length)
{
  while(length > 0)
    {
      *dest = *(src+position);
      dest++;
      src++;
      length--;
    }
} 

char codon_table(char *temp){

  char pep;
  if (strcmp(temp,"GCT")==0||strcmp(temp,"GCC")==0||strcmp(temp,"GCA")==0||strcmp(temp,"GCG")==0){
    pep = 'A';
  }else if (strcmp(temp,"CGT")==0||strcmp(temp,"CGC")==0||strcmp(temp,"CGA")==0||strcmp(temp,"CGG")==0||strcmp(temp,"AGA")==0||strcmp(temp,"AGG")==0){ 
    pep = 'R';
  }else if (strcmp(temp,"AAT")==0||strcmp(temp,"AAC")==0){
    pep = 'N';
  }else if (strcmp(temp,"GAT")==0||strcmp(temp,"GAC")==0){
    pep = 'D';
  }else if (strcmp(temp,"TGT")==0||strcmp(temp,"TGC")==0){
    pep = 'C';
  }else if (strcmp(temp,"CAA")==0||strcmp(temp,"CAG")==0){
    pep = 'Q';
  }else if (strcmp(temp,"GAA")==0||strcmp(temp,"GAG")==0){
    pep = 'E';
  }else if (strcmp(temp,"GGT")==0||strcmp(temp,"GGC")==0||strcmp(temp,"GGA")==0||strcmp(temp,"GGG")==0){
    pep = 'G';
  }else if (strcmp(temp,"CAT")==0||strcmp(temp,"CAC")==0){
    pep = 'H';
  }else if (strcmp(temp,"ATT")==0||strcmp(temp,"ATC")==0||strcmp(temp,"ATA")==0){
    pep = 'I';
  }else if (strcmp(temp,"ATG")==0){
    pep = 'M';
  }else if (strcmp(temp,"TTA")==0||strcmp(temp,"TTG")==0||strcmp(temp,"CTT")==0||strcmp(temp,"CTC")==0||strcmp(temp,"CTA")==0||strcmp(temp,"CTG")==0){ 
    pep = 'L';
  }else if (strcmp(temp,"AAA")==0||strcmp(temp,"AAG")==0){
    pep = 'K';
  }else if (strcmp(temp,"ATG")==0){
    pep = 'M';
  }else if (strcmp(temp,"TTT")==0||strcmp(temp,"TTC")==0){
    pep = 'F';
  }else if (strcmp(temp,"CCT")==0||strcmp(temp,"CCC")==0||strcmp(temp,"CCA")==0||strcmp(temp,"CCG")==0){
    pep = 'P';
  }else if (strcmp(temp,"TCT")==0||strcmp(temp,"TCC")==0||strcmp(temp,"TCA")==0||strcmp(temp,"TCG")==0||strcmp(temp,"AGT")==0||strcmp(temp,"AGC")==0){ 
    pep = 'S';
  }else if (strcmp(temp,"ACT")==0||strcmp(temp,"ACC")==0||strcmp(temp,"ACA")==0||strcmp(temp,"ACG")==0){
    pep = 'T';
  }else if (strcmp(temp,"TGG")==0){
    pep = 'W';
  }else if (strcmp(temp,"TAT")==0||strcmp(temp,"TAC")==0){
    pep = 'Y';
  }else if (strcmp(temp,"GTT")==0||strcmp(temp,"GTC")==0||strcmp(temp,"GTA")==0||strcmp(temp,"GTG")==0){
    pep = 'V';
  }else if (strcmp(temp,"TAG")==0||strcmp(temp,"TGA")==0||strcmp(temp,"TAA")==0){  
    pep = '*';
  }else {
    if (temp[0] == 'N' || temp[1] == 'N' || temp[2] == 'N'){
      pep = 'X';
    }else{
      pep = 'X';
      /* printf("nonsense character %s in fasta file\n",temp); */
    }
  }
  return pep;
}


void translate_seq(char *O, FILE *fp, char *head){

    char *temp;
    int i,j;
    int count = 0;
    temp = (char *)malloc(3 * sizeof(char));
    int seq_len = strlen(O);
   

    fprintf(fp, ">%s_1\n", head);
    for (i=0; i<seq_len-3; i+=3){
      mysubstr(temp, O, i, 3);
      for(j=0; j<3; j++){
	if ((int)temp[j] > 90){
	  temp[j] -= 32;
	}
      }

      fprintf(fp,"%c",codon_table(temp));
      count ++;
      if (count % 60 == 0 ){
	fprintf(fp, "\n");
      } 
    } 
    count=0;
    fprintf(fp, "\n>%s_2\n", head);
    for (i=1; i<seq_len-3; i+=3){
      mysubstr(temp, O, i, 3);
      for(j=0; j<3; j++){
	if ((int)temp[j] > 90){
	  temp[j] -= 32;
	}
      }

      fprintf(fp,"%c",codon_table(temp));
      count ++;
      if (count % 60 == 0 ){
	fprintf(fp, "\n");
      } 
    } 
    count=0;
    fprintf(fp, "\n>%s_3\n", head);
    for (i=2; i<seq_len-3; i+=3){
      mysubstr(temp, O, i, 3);
      for(j=0; j<3; j++){
	if ((int)temp[j] > 90){
	  temp[j] -= 32;
	}
      }

      fprintf(fp,"%c",codon_table(temp));
      count ++;
      if (count % 60 == 0 ){
	fprintf(fp, "\n");
      } 
    } 


}

int main (int argc, char **argv){

  FILE *fp;
  int num_seqs;
  char **obs_seqs = NULL;
  char *seq_file = NULL;
  char *pep_file = NULL;
  char *head = NULL;
  char c;

  /* read command line argument */
  while ((c=getopt(argc, argv, "d:h:p:")) != -1){

    switch (c){
    case 'd':
      seq_file = optarg;   /* DNA file */
      break;
      
    case 'h':
      head = optarg;   /* DNA SEQ name */
      break;      
 
    case 'p':
      pep_file = optarg;   /* PEP file */
      break;      
   
    }
  }

 /* read observation seq */
  fp = fopen(seq_file, "r");
  read_obs_seq_from_file(fp, &num_seqs, &obs_seqs);
  fclose(fp);

  /* translate */
  fp = fopen(pep_file, "w");
  translate_seq(obs_seqs[0], fp, head);
  fclose(fp);

  return 0;
}

