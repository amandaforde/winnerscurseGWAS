#!/bin/sh
#SBATCH -J T2DB
#SBATCH --partition=highmem

for i in {1..22}; do

sbatch -p highmem -N 1 -n 16 T2DB.sh $i;

done
