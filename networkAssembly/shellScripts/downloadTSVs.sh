#! /bin/bash

for TF in `cat ./analysisList.tab | grep hg19 | cut -f 1`; do file="http://dbarchive.biosciencedbc.jp/kyushu-u/hg19/target/"$TF".5.tsv"; touch $TF; curl $file > $TF; done