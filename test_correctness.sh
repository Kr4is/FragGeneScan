#!/bin/bash

echo --------- Clean and Make ---------
echo ""

echo "---> Make clean"
make clean
echo ""

echo "---> Make fgs"
make fgs
echo ==================================================
echo ""

echo --------- Clean Test Directory ---------
echo ""

echo "---> Remove tests directory"
rm -r ./tests
echo ""

echo "---> Create tests directory"
mkdir tests
echo ==================================================
echo ""

echo --------- Executions ---------
echo ""

echo "---> Complete genome secuence"
./run_FragGeneScan.pl -genome=./example/NC_000913.fna -out=./tests/NC_000913-fgs  -complete=1  -train=complete
echo ""

echo "---> Sequencing reads"
./run_FragGeneScan.pl -genome=./example/NC_000913-454.fna -out=./tests/NC_000913-454-fgs  -complete=0  -train=454_10
echo ""

echo "---> Assembly contigs"
./run_FragGeneScan.pl -genome=./example/contigs.fna -out=./tests/contigs-fgs  -complete=1  -train=complete
echo ==================================================
echo ""

echo --------- Test ---------
echo ""

echo "---> Complete genome secuence"
diff ./example/NC_000913-fgs.out ./tests/NC_000913-fgs.out
echo ""

echo "---> Sequencing reads"
diff ./example/NC_000913-454-fgs.out ./tests/NC_000913-454-fgs.out
echo ""

echo "---> Assembly contigs"
diff ./example/contigs-fgs.out ./tests/contigs-fgs.out
echo ==================================================