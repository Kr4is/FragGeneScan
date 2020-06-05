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

echo --------- Clean Project ---------
echo ""

echo "---> Make clean"
make clean
echo ==================================================

echo --------- Make Project ---------
echo ""

echo "---> Make fgs"
make CC=mpicc fgs
echo ==================================================
echo ""