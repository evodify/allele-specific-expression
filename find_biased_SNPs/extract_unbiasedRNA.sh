#!/bin/sh

sed 's/"//g;s/ /\t/g' $1 | cut -f 2 > $1_Names.csv
grep -Fwf $1_Names.csv $2 > $2_unbiased.csv
rm $1_Names.csv
