#!/bin/bash
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
if [[ $# -eq 0 ]] ; then
    echo 'Set the number of cores to test'
    echo 'Example: test_correctess_mpi.sh num_cores'
    exit 1
fi

echo --------- Make Project ---------
echo ""

echo "---> Make fgs"
make CC=mpicc fgs
echo ==================================================
echo ""

echo --------- Test Directory Creation ---------
echo ""

echo "---> Create tests directory"
mkdir tests
echo "tests directory created"
echo ==================================================
echo ""

echo --------- Executions ---------
echo ""

echo "---> Complete genome secuence"
./run_FragGeneScan_mpi.pl -genome=./example/NC_000913.fna -out=./tests/NC_000913-fgs  -complete=1  -train=complete -thread=$1
echo ""

echo "---> Sequencing reads"
./run_FragGeneScan_mpi.pl -genome=./example/NC_000913-454.fna -out=./tests/NC_000913-454-fgs  -complete=0  -train=454_10 -thread=$1
echo ""

echo "---> Assembly contigs"
./run_FragGeneScan_mpi.pl -genome=./example/contigs.fna -out=./tests/contigs-fgs  -complete=1  -train=complete -thread=$1
echo ==================================================
echo ""

echo --------- Test ---------
echo ""

echo "---> Complete genome secuence"
DIFF=$(diff ./example/NC_000913-fgs.out ./tests/NC_000913-fgs.out)
if [ "$DIFF" == "" ] 
then
    echo "Test Passed"
else
    echo "TEST ERROR"
fi
echo ""

echo "---> Sequencing reads"
DIFF=$(diff ./example/NC_000913-454-fgs.out ./tests/NC_000913-454-fgs.out)
if [ "$DIFF" == "" ] 
then
    echo "Test Passed"
else
    echo "TEST ERROR"
fi
echo ""

echo "---> Assembly contigs"
DIFF=$(diff ./example/contigs-fgs.out ./tests/contigs-fgs.out)
if [ "$DIFF" == "" ] 
then
    echo "Test Passed"
else
    echo "TEST ERROR"
fi
echo ==================================================
echo ""

echo --------- Test Directory Suppression ---------
echo ""

echo "---> Remove tests directory"
rm -r ./tests
echo ==================================================
echo ""

echo --------- Clean Project ---------
echo ""

echo "---> Make clean"
make clean
echo ==================================================