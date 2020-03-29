/*
  FragGeneScan: predicting genes in short and error-prone reads.
	Copyright © 2010 Mina Rho, Yuzhen Ye and Haixu Tang.
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

double **dmatrix(int num_row, int num_col);
double *dvector(int nh);
int **imatrix(int num_row, int num_col);
int *ivector(int nh);

void free_dvector(double *v);
void free_dmatrix(double **m,int num_row);
void free_ivector(int *v);
void free_imatrix(int **m,int num_row);

int tr2int (char *nt);
int nt2int (char nt);
int nt2int_rc (char nt);


int trinucleotide (char a, char b, char c);
void get_rc_dna_indel(char *dna, char *dna1);
double log2(double a);
void get_protein(char *dna, char *protein, int strand, int whole_genome);
void print_usage();
