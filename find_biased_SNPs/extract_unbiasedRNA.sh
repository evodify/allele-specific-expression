#!/bin/sh

# to extract the RNA data SNPs that showed no bias in the DNA data.
#
# command:
# 
# sh extractUnbiased.sh DNAinput1_unbiased.csv RNAinput.csv

sed 's/"//g;s/ /\t/g' $1 | cut -f 2 > $1_Names.csv
grep -Fwf $1_Names.csv $2 > $2_unbiased.csv
rm $1_Names.csv
