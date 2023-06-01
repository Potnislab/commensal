#!/bin/bash


#!/bin/bash

# (C) Copyright 2017 Isai Salas Gonzalez
#
#    This file is part of gfobap.
#
#    gfobap is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    gfobap is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with gfobap.  If not, see <http://www.gnu.org/licenses/>.

#This wrapper should be run in the folder where the input matrix is located
#Mofidy the spath  to wherever the enrichment scripts are located
spath=/scratch/aubnxp/biosample/workflow/enrichment_tests;
file=$1;

module load R/4.1.0
#Run phyloglm
echo "Running phyloglm";
Rscript --vanilla "$spath"/oh.phyloglm.R splitted_file_1.tsv
