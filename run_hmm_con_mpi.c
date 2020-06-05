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

#include "mpi.h"

#define STRINGLEN 4096

int main (int argc, char **argv)
{
  clock_t start = clock();
  int i, j, c;
  HMM hmm;
  TRAIN train;
  int wholegenome;
  int format = 0;
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
  int count = 0;
  int total = 0;
  char mystring[STRINGLEN+4] = "";
  int *obs_seq_len;
  int *seq_start_pointers;
  int seq_pointer;
  int bp_count;  /* count the length of each line in input file */

  // Viterbi variables
  FILE *out;
	FILE *aa;
	FILE *dna;
	char *obs_head;
	char *obs_seq;
	int cg;

  int threadnum = 1;
  int numprocs, myid;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  if(myid == 0){
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
        if (access(seq_file, F_OK) == -1){
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
        if (access(hmm_file, F_OK) == -1){
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
    if (access(mstate_file, F_OK) == -1){
      fprintf(stderr, "Forward prob. file [%s] does not exist\n", mstate_file);
      exit(1);
    }
    if (access(rstate_file, F_OK) == -1){
      fprintf(stderr, "Backward prob. file [%s] does not exist\n", rstate_file);
      exit(1);
    }
    if (access(nstate_file, F_OK) == -1){
      fprintf(stderr, "noncoding prob. file [%s] does not exist\n", nstate_file);
      exit(1);
    }
    if (access(sstate_file, F_OK) == -1){
      fprintf(stderr, "start prob. file [%s] does not exist\n", sstate_file);
      exit(1);
    }
    if (access(pstate_file, F_OK) == -1){
      fprintf(stderr, "stop prob. file [%s] does not exist\n", pstate_file);
      exit(1);
    }
    if (access(s1state_file, F_OK) == -1){
      fprintf(stderr, "start1 prob. file [%s] does not exist\n", s1state_file);
      exit(1);
    }
    if (access(p1state_file, F_OK) == -1){
      fprintf(stderr, "stop1 prob. file [%s] does not exist\n", p1state_file);
      exit(1);
    }
    if (access(dstate_file, F_OK) == -1){
      fprintf(stderr, "pwm dist. file [%s] does not exist\n", dstate_file);
      exit(1);
    }
    if (access(hmm_file, F_OK) == -1){
      fprintf(stderr, "hmm file [%s] does not exist\n", hmm_file);
      exit(1);
    }
    
    /* read all initial model */
    hmm.N=NUM_STATE;
    get_train_from_file(hmm_file, &hmm, mstate_file, rstate_file, nstate_file, sstate_file, pstate_file, s1state_file, p1state_file, dstate_file, &train);
  }

  MPI_Bcast(&hmm, sizeof(HMM), MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(&train, sizeof(TRAIN), MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(wholegenome, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(format, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(out_header, strlen(out_header), MPI_CHAR, 0, MPI_COMM_WORLD);
  MPI_Bcast(seq_file, strlen(seq_file), MPI_CHAR, 0, MPI_COMM_WORLD);

  // Initialize thread data structure
  sprintf(mystring, "%s.out.tmp%d", out_header, myid);
  out = fopen(mystring, "w");
  sprintf(mystring, "%s.faa.tmp%d", out_header, myid);
  aa = fopen(mystring, "w");
  sprintf(mystring, "%s.ffn.tmp%d", out_header, myid);
  dna = fopen(mystring, "w");

  if (myid == 0){
    fp = fopen(seq_file, "r");
    while (fgets (mystring, sizeof(mystring), fp)){
      if (mystring[0] == '>'){
        count++;
      }
    }
    printf("no. of seqs: %d\n", count);
    // rewind file to de top
    rewind(fp);
  }
  MPI_Bcast(count, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Sequences length
  obs_seq_len = (int *) malloc(count * sizeof(int));
  // Sequences start pointers
  seq_start_pointers = (int *) malloc(count * sizeof(int));

  if (myid == 0){
    seq_pointer = 0;
    i = 0;
    count = 0;
    while (fgets(mystring, sizeof(mystring), fp)){
      if (mystring[0] == '>'){
        if (i > 0){
          // previous sequence element
          obs_seq_len[count] = i;
          count++;
        }
        seq_start_pointers[count] = seq_pointer;
        i = 0;
      }else{
        bp_count = strlen(mystring);
        // chr(13) => "\r"
        // chr(10) => "\n"
        while(mystring[bp_count - 1] == 10 || mystring[bp_count - 1] == 13){
          bp_count--;
        }
        i += bp_count;
      }
      seq_pointer = ftell(fp);
    }
    // Last element
    obs_seq_len[count] = i;
    // rewind file to de top
    rewind(fp);
  }
  MPI_Bcast(obs_seq_len, count, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(seq_start_pointers, count, MPI_INT, 0, MPI_COMM_WORLD);

  int my_start;
  for(i=0; i<myid; i++){
    if(i < (count % numprocs)){
      my_start += count + 1;
    }else{
      my_start += count;
    }
  }
  int mi_size = myid < (count % numprocs) ? (count/numprocs) + 1 : count/numprocs;

  //count % numprocs;
  //count / numprocs;
  total = 0;
  count = 0;
  j = 0;

  while (!(feof(fp))){
    memset(mystring, '\0', sizeof(mystring));
    fgets (mystring, sizeof(mystring), fp);
    bp_count = strlen(mystring);
    // chr(13) => "\r"
    // chr(10) => "\n"
    while(mystring[bp_count - 1] == 10 || mystring[bp_count - 1] == 13){
      bp_count--;
    }

    if (mystring[0] == '>' || feof(fp)){
      if (feof(fp)){
        memcpy(obs_seq + j, mystring, bp_count);
        j += bp_count;
      }
      if ((count > 0) || feof(fp)){
        // Deal with the thread
        cg = get_prob_from_cg(&hmm, &train, obs_seq); // cg - 26 Ye April 16, 2016
        if (strlen(obs_seq)>70){
          viterbi(&hmm, &train, obs_seq, out, aa, dna, obs_head, wholegenome, cg, format);
        }
        free(obs_head);
        free(obs_seq);
        obs_head = NULL;
        obs_seq = NULL;
        count = 0;
      }

      if (!(feof(fp))){
        obs_head = (char *) malloc((bp_count + 1) * sizeof(char));
        memset(obs_head, 0, (bp_count + 1) * sizeof(char));
        memcpy(obs_head, mystring, bp_count);

        obs_seq = (char*) malloc((obs_seq_len[total] + 1) * sizeof(char));
        memset(obs_seq, '\0', (obs_seq_len[total] + 1) * sizeof(char));
        total++;
        count++;
        j = 0;
      }

    }else{
      memcpy(obs_seq + j, mystring, bp_count);
      j += bp_count;
    }
    if (feof(fp)){
      break;
    }
  }
  
  fclose(out);
  fclose(aa);
  fclose(dna);

  MPI_Finalize();

  clock_t end = clock();
  printf("Clock time used (by %d proceses = %.2f mins\n", numprocs, (end - start) / (60.0 * CLOCKS_PER_SEC));
}

