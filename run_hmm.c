/*
  FragGeneScan: predicting genes in short and error-prone reads.
	Copyright © 2010-2018 Mina Rho, Yuzhen Ye and Haixu Tang.
  Copyright © 2020 Bruno Cabado Lousa.
	
	This file is part of FragGeneScan.

  FragGeneScan is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FragGeneScan is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with FragGeneScan.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>
#include <time.h>
#include "hmm.h"
#include "util_lib.h"

#define STRINGLEN 4096

typedef struct thread_data
{
	FILE *out;
	FILE *aa;
	FILE *dna;
	char *obs_head;
	char *obs_seq;
	int wholegenome;
	int cg;
	int format;
	HMM *hmm;
	TRAIN *train;
} thread_data;

int main (int argc, char **argv)
{
  clock_t start = clock();
  int i, j, c;
  HMM hmm;
  TRAIN train;
  int wholegenome;
  int format=0;
  FILE *fp;
  char hmm_file[STRINGLEN] = "";
  char out_header[STRINGLEN] = "";
  char seq_file[STRINGLEN] = "";
  char train_file[STRINGLEN] = "";
  char mstate_file[STRINGLEN] = "";
  char rstate_file[STRINGLEN] = "";
  char nstate_file[STRINGLEN] = "";
  char sstate_file[STRINGLEN] = "";
  char pstate_file[STRINGLEN] = "";
  char s1state_file[STRINGLEN] = "";     /* stop codon of gene in - stand */
  char p1state_file[STRINGLEN] = "";
  char dstate_file[STRINGLEN] = "";
  char train_dir[STRINGLEN] = "";
  int count=0;
  int total = 0;
  char mystring[STRINGLEN+4] = "";
  int *obs_seq_len;
  int bp_count;  /* count the length of each line in input file */

  int threadnum = 1;

  strncpy(train_dir, argv[0], strlen(argv[0])-12);
  strcat(train_dir, "train/");
  strcpy(mstate_file, train_dir);
  strcat(mstate_file, "gene");
  strcpy(rstate_file, train_dir);
  strcat(rstate_file, "rgene");
  strcpy(nstate_file, train_dir);
  strcat(nstate_file, "noncoding");
  strcpy(sstate_file, train_dir);
  strcat(sstate_file, "start");
  strcpy(pstate_file, train_dir);
  strcat(pstate_file, "stop");
  strcpy(s1state_file, train_dir);
  strcat(s1state_file, "stop1");
  strcpy(p1state_file, train_dir);
  strcat(p1state_file, "start1");
  strcpy(dstate_file, train_dir);
  strcat(dstate_file, "pwm");


  /* read command line argument */
  if (argc <= 8){    
    fprintf(stderr, "ERROR: You missed some parameters for input\n");
    print_usage();
    exit(EXIT_FAILURE);
  }

  while ((c=getopt(argc, argv, "fs:o:w:t:p:")) != -1){
    switch (c){
    case 's':
      strcpy(seq_file, optarg);
      if (access(seq_file, F_OK)==-1){
        fprintf(stderr, "ERROR: Sequence file [%s] does not exist\n", seq_file);
        print_usage();
        exit(EXIT_FAILURE);
      }
      break;  
    case 'w':
      wholegenome = atoi(optarg);
      if (wholegenome != 0 && wholegenome != 1){
        fprintf(stderr, "ERROR: An incorrect value for the option -w was entered\n");
        print_usage();
        exit(EXIT_FAILURE);
      }
      break;
    case 'p':
      threadnum = atoi(optarg);
      if (threadnum < 1){
        fprintf(stderr, "ERROR: An incorrect value [%d] for the option -p was entered\n", threadnum);
        print_usage();
        exit(EXIT_FAILURE);
      }
      printf("Using %d threads.\n", threadnum);
      break;
    case 'o':
      strcpy(out_header, optarg);
      break;
    case 't':
      strcpy(train_file, optarg);
      strcpy(hmm_file, train_dir);
      strcat(hmm_file, train_file);

      if (access(hmm_file, F_OK)==-1){
        fprintf(stderr, "ERROR: The file for model parameters [%s] does not exist\n", hmm_file);
        print_usage();
        exit(EXIT_FAILURE);
      }
      break;
    case 'f':
      format = 1;
      break;
    }
  }

  
  /* check whether the specified files exist */
  if (access(mstate_file, F_OK)==-1){
    fprintf(stderr, "Forward prob. file [%s] does not exist\n", mstate_file);
    exit(1);
  }
  if (access(rstate_file, F_OK)==-1){
    fprintf(stderr, "Backward prob. file [%s] does not exist\n", rstate_file);
    exit(1);
  }
  if (access(nstate_file, F_OK)==-1){
    fprintf(stderr, "noncoding prob. file [%s] does not exist\n", nstate_file);
    exit(1);
  }
  if (access(sstate_file, F_OK)==-1){
    fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
    exit(1);
  }
  if (access(pstate_file, F_OK)==-1){
    fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
    exit(1);
  }
  if (access(s1state_file, F_OK)==-1){
    fprintf(stderr, "start1 prob. file [%s] does not exist\n", s1state_file);
    exit(1);
  }
  if (access(p1state_file, F_OK)==-1){
    fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
    exit(1);
  }
  if (access(dstate_file, F_OK)==-1){
    fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
    exit(1);
  }
  if (access(hmm_file, F_OK)==-1){
    fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
    exit(1);
  }
  
  /* read all initial model */
  hmm.N=NUM_STATE;
  get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file,s1state_file, p1state_file, dstate_file, &train);

  // Initialize thread data structure
  threadnum = 1;
  thread_data data;

  sprintf(mystring, "%s.out", out_header); 
  data.out = fopen(mystring, "w");
  sprintf(mystring, "%s.faa", out_header);
  data.aa = fopen(mystring, "w");
  sprintf(mystring, "%s.ffn", out_header);
  data.dna = fopen(mystring, "w");

  data.hmm = (HMM*)malloc(sizeof(HMM));
  memcpy(data.hmm, &hmm, sizeof(HMM));
  data.train = (TRAIN*)malloc(sizeof(TRAIN));
  memcpy(data.train, &train, sizeof(TRAIN));

  data.wholegenome = wholegenome;
  data.format = format;

  void *status;
  fp = fopen (seq_file, "r");
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
      count++;
    }
  }
  obs_seq_len = (int *)malloc(count * sizeof(int));
  printf("no. of seqs: %d\n", count);  

  i = 0;
  count = 0;
  rewind(fp);
  while ( fgets (mystring , sizeof mystring , fp) ){
    if (mystring[0] == '>'){
      if (i>0){
        obs_seq_len[count] = i;
        count++;
      }
      i = 0;
    }else{
      bp_count = strlen(mystring);
      while(mystring[bp_count-1] == 10 || mystring[bp_count-1]==13){
	      bp_count --;
      }
      i += bp_count;
    }
  }
  obs_seq_len[count] = i;

  rewind(fp);
  total = 0;
  count = 0;
  j = 0;

  while (!(feof(fp)))
  {
    memset(mystring, '\0', sizeof mystring);
    fgets (mystring , sizeof mystring  , fp);
    bp_count = strlen(mystring);
    while(mystring[bp_count - 1] == 10 || mystring[bp_count - 1]==13){
      bp_count --;
    }

    if (mystring[0] == '>' || feof(fp)){
      if (feof(fp))
      {
        memcpy(data.obs_seq + j, mystring, bp_count);
        j += bp_count;
      }
      if ((count > 0 && count % threadnum == 0) || feof(fp))
      {
        // Deal with the thread
        data.cg = get_prob_from_cg(data.hmm, data.train, data.obs_seq); //cg - 26 Ye April 16, 2016
        if (strlen(data.obs_seq)>70){
          viterbi(data.hmm, data.train, data.obs_seq, data.out, data.aa, data.dna, data.obs_head, data.wholegenome, data.cg, data.format);
        }

        free(data.obs_head);
        free(data.obs_seq);
        data.obs_head = NULL;
        data.obs_seq = NULL;
        count = 0;
      }

      if (!(feof(fp)))
      {
        data.obs_head = (char *)malloc((bp_count+1) * sizeof(char));
        memset(data.obs_head, 0, (bp_count+1) * sizeof(char));
        memcpy(data.obs_head, mystring, bp_count);

        data.obs_seq = (char*)malloc((obs_seq_len[total] + 1) * sizeof(char));
        memset(data.obs_seq, '\0', (obs_seq_len[total] + 1) * sizeof(char));
        total++;
        count++;
        j = 0;
      }

    }else{
      memcpy(data.obs_seq + j, mystring, bp_count);
      j += bp_count;
    }
    if (feof(fp))
    {
      break;
    }
  }
  fclose(data.out);
  fclose(data.aa);
  fclose(data.dna);

  clock_t end = clock();
  printf("Clock time used (by %d threads) = %.2f mins\n", threadnum, (end - start) / (60.0 * CLOCKS_PER_SEC));
}

