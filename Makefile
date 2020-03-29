#
#   FragGeneScan: predicting genes in short and error-prone reads.
#   Copyright © 2020 Bruno Cabado Lousa.
#	
#   This file is part of FragGeneScan.
#
#   FragGeneScan is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   FragGeneScan is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with FragGeneScan.  If not, see <https://www.gnu.org/licenses/>.
#

CC=	gcc
CFLAG= -O3
SRCS=	util_lib.c hmm_lib.c run_hmm.c
OBJ=	util_lib.o hmm_lib.o run_hmm.o 

hmm.obj:	$(SRCS)
	$(CC) $(CFLAG) -c $(SRCS)

fgs:	$(OBJ)
	$(CC)  $(CFLAG) -o FragGeneScan util_lib.o hmm_lib.o run_hmm.o  -lm

clean:
	rm -rf *.o FragGeneScan* *~
