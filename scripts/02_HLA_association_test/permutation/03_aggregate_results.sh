#! /bin/bash

base_dir=""
cd ${base_dir}

echo -e "P\tFN" > permutation/results/permutation_summary.TCR.1_200.result.txt

for i in {1..200}
do
    awk 'BEGIN {OFS="\t"}{if(FNR!=1)print $6,FILENAME}' permutation/results/TCR.${i}/permutation.TCR.${i}.*.result.txt >> permutation/results/permutation_summary.TCR.1_200.result.txt
done
gzip permutation/results/permutation_summary.TCR.1_200.result.txt
