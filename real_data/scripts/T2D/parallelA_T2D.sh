#!/bin/sh
#SBATCH -J T2DA
#SBATCH --partition=highmem

for i in {1..22}; do

sbatch -p highmem -N 1 -n 16 T2DA.sh $i;

done
