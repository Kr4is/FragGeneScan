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