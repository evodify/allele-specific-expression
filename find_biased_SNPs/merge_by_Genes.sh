#!/bin/sh

# To merge results of ASE analyses (either `*_signASE.csv` or `*_allASE.csv`) by gene names.
#
# command:
#
# sh merge_ASE_by_Gene.sh RNAinput1_unbiased_allASE.csv

# ratio per SNP averaged per gene:
awk  'BEGIN {print "Gene\tnumber_of_SNPs\tproportion_of_homeologueA"}; NR>1, NF>1 {a[$1] += $3/($3+$4); n[$1] = n[$1]+=1}; END {for (i in a) {print i"\t"n[i]"\t"(a[i]/n[i])} }' $1 > $1.Genes

# ratio per gene:
# awk 'BEGIN {print "Gene\tnumber_of_SNPs\thomeologueA_count\thomeologueB_count\tproportion_of_homeologueA"}; NR>1, NF>1 {a[$1] = a[$1]+=$3; b[$1] = b[$1]+=$4; n[$1] = n[$1]+=1}; END {for (i in a) {print i"\t"n[i]"\t"a[i]"\t"b[i]"\t"(a[i]/(a[i]+b[i]))} }' $1 > $1.Genes
