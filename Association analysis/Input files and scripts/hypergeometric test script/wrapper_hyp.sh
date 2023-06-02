#!/bin/bash


#This wrapper should be run in the folder where the input matrix is located
#Mofidy the spath  to wherever the enrichment scripts are located
spath=/scratch/aubnxp/biosample/enrichment_tests;
file=$1;

module load R/4.1.0

#Run the hypergeometric test
echo "Running hypergeometric";
Rscript --vanilla "$spath"/hypergeometric.R matrix_raw_orthogroups.tsv
Rscript --vanilla "$spath"/hypergeometric.R matrix_binary_orthogroups.tsv