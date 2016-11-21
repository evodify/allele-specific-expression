#!/bin/sh

# To merge results of ASE analyses (either `*_signASE.csv` or `*_allASE.csv`)
# with all the expression information (`Genes_*` file)
#
# command:
#
# sh merge_ASE_Gene.sh RNAinput1_unbiased_allASE.csv Genes_RNAinput1.csv

sort $1 > $1_temp
sort $2 > $2_temp
join -j 1 $1_temp $2_temp > $2_mergedASE.csv
rm *_temp
