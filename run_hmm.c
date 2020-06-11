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
  int num_sequences = 0;
  int count = 0;
  int current_seq = 0;
  char mystring[STRINGLEN] = "";
  char filename[STRINGLEN+16] = "";
  int *obs_seq_len, *start_sequences, *end_sequences;
  long int *seq_start_pointers;
  long int seq_pointer, num_start_seq, num_final_seq;
  int next_proc;
  int bp_count;  /* count the length of each line in input file */

  // Viterbi variables
  FILE *out;
	FILE *aa;
	FILE *dna;
	char *obs_head;
	char *obs_seq;
	int cg;

  // MPI variables
  int myid, num_procs;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

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

  while ((c=getopt(argc, argv, "fs:o:w:t:")) != -1){
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

  // Create output temp files for each process
  sprintf(filename, "%s.out.%d", out_header, myid);
  out = fopen(filename, "w");
  sprintf(filename, "%s.faa.%d", out_header, myid);
  aa = fopen(filename, "w");
  sprintf(filename, "%s.ffn.%d", out_header, myid);
  dna = fopen(filename, "w");

  // Open sequences file
  fp = fopen(seq_file, "r");

  if(myid == 0){
    // Get the number of sequences in sequences file
    while (fgets(mystring, sizeof(mystring), fp)){
      if (mystring[0] == '>'){
        num_sequences++;
      }
    }
    printf("no. of seqs: %d\n", num_sequences);

    // Rewind file to de top
    rewind(fp);
  }
  /* Broadcast sequences number from process 0 to all processes */
  MPI_Bcast(&num_sequences, 1, MPI_INT, 0, MPI_COMM_WORLD);

  // Sequences lengths
  obs_seq_len = (int *) malloc(num_sequences * sizeof(int));
  // Sequences start positions in sequences file
  seq_start_pointers = (long int *) malloc(num_sequences * sizeof(long int));
  // Start sequences for each process
  start_sequences = (int *) malloc(num_procs * sizeof(int));
  // End sequences for each process
  end_sequences = (int *) malloc(num_procs * sizeof(int));

  long int total_seqs_len = 0;

  if(myid == 0){
    seq_pointer = 0;
    i = 0;
    count = 0;
    while (fgets(mystring, sizeof(mystring), fp)){
      if (mystring[0] == '>'){
        if (i > 0){
          // previous sequence element
          obs_seq_len[count] = i;
          total_seqs_len += i;
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
    total_seqs_len += i;
    // Rewind file to de top
    rewind(fp);
    
    int seqs_len_per_proc = total_seqs_len / num_procs;
    int curr_proc = 0;
    long int length_accum = 0;

    for(current_seq = 0; current_seq < num_sequences; current_seq++){
      if(length_accum >= seqs_len_per_proc * curr_proc){
        if(curr_proc == 0){
          end_sequences[num_procs - 1] = num_sequences;
        }else{
          end_sequences[curr_proc - 1] = current_seq;
        }
        start_sequences[curr_proc] = current_seq;
        curr_proc++;
      }
      length_accum += obs_seq_len[current_seq];
    }
  }
  
  /* Broadcast sequences lengths from process 0 to all processes */
  MPI_Bcast(obs_seq_len, num_sequences, MPI_INT, 0, MPI_COMM_WORLD);
  /* Broadcast sequences start positions in sequences file from process 0 to all processes */
  MPI_Bcast(seq_start_pointers, num_sequences, MPI_LONG, 0, MPI_COMM_WORLD);
  /* Broadcast start sequences from process 0 to all processes */
  MPI_Bcast(start_sequences, num_procs, MPI_INT, 0, MPI_COMM_WORLD);
  /* Broadcast end sequences from process 0 to all processes */
  MPI_Bcast(end_sequences, num_procs, MPI_INT, 0, MPI_COMM_WORLD);

  // Get sequences to do by each process
  num_start_seq = start_sequences[myid];
  num_final_seq = end_sequences[myid];

  current_seq = num_start_seq;
  count = 0;
  j = 0;

  // Move process to each initial position in file
  fseek(fp, seq_start_pointers[num_start_seq], SEEK_SET);

  // Calculate process sequences
  while (!(feof(fp)) && (current_seq <= num_final_seq)){
    memset(mystring, '\0', sizeof(mystring));
    fgets(mystring, sizeof(mystring), fp);
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

        obs_seq = (char*) malloc((obs_seq_len[current_seq] + 1) * sizeof(char));
        memset(obs_seq, '\0', (obs_seq_len[current_seq] + 1) * sizeof(char));
        current_seq++;
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
  
  // Close files
  fclose(out);
  fclose(aa);
  fclose(dna);

  // Wait for all processes computation
  MPI_Barrier(MPI_COMM_WORLD);

  if(myid == 0){
    // Create final output files
    FILE *current_file;
    int current_proc;
    // out file
    FILE *final_out;
    sprintf(filename, "%s.out", out_header);
    final_out = fopen(filename, "w");
    for(current_proc = 0; current_proc < num_procs; current_proc++){
      sprintf(filename, "%s.out.%d", out_header, current_proc);
      current_file = fopen(filename, "r");
      while(fgets(mystring, sizeof(mystring), current_file)) {
        fprintf(final_out, "%s", mystring);
      }
      fclose(current_file);
      // Remove temp file
      remove(filename);
    }
    fclose(final_out);
    // aa file
    FILE *final_aa;
    sprintf(filename, "%s.faa", out_header);
    final_aa = fopen(filename, "w");
    for(current_proc = 0; current_proc < num_procs; current_proc++){
      sprintf(filename, "%s.faa.%d", out_header, current_proc);
      current_file = fopen(filename, "r");
      while(fgets(mystring, sizeof(mystring), current_file)) {
        fprintf(final_aa, "%s", mystring);
      }
      fclose(current_file);
      // Remove temp file
      remove(filename);
    }
    fclose(final_aa);
    // dna file
    FILE *final_dna;
    sprintf(filename, "%s.ffn", out_header);
    final_dna = fopen(filename, "w");
    for(current_proc = 0; current_proc < num_procs; current_proc++){
      sprintf(filename, "%s.ffn.%d", out_header, current_proc);
      current_file = fopen(filename, "r");
      while(fgets(mystring, sizeof(mystring), current_file)) {
        fprintf(final_dna, "%s", mystring);
      }
      fclose(current_file);
      // Remove temp file
      remove(filename);
    }
    fclose(final_dna);
  }

  clock_t end = clock();
  if(myid == 0){
    printf("Clock time used (by %d processes) = %.2f mins\n", num_procs, (end - start) / (60.0 * CLOCKS_PER_SEC));
  }
  MPI_Finalize();
}

