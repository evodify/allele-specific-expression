#!/bin/sh

# To merge results of ASE analyses (either `*_signASE.csv` or `*_allASE.csv`)
# with all the expression information (`Genes_*` file)
#
# command:
#
# sh merge_ASE_with_Gene.sh RNAinput1_unbiased_allASE.csv Genes_RNAinput1.csv

echo -e "" > header_temp
sort $1 > $1_temp
sort $2 > $2_temp
join -j 1 $1_temp $2_temp | sed 's/ /\t/g;1i\gene\tASE_prob\tnumber_of_SNPs\thomeologueA_count\thomeologueB_count\tproportion_of_homeologueA' > $2_mergedASE.csv
rm *_temp
