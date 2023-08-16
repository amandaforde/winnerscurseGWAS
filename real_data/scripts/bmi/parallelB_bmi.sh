#!/bin/sh
#SBATCH -J bmiB
#SBATCH --partition=highmem

for i in {1..22}; do

sbatch -p highmem -N 1 -n 16 bmiB.sh $i;

done
