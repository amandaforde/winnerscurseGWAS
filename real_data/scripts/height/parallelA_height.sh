#!/bin/sh
#SBATCH -J heightA
#SBATCH --partition=highmem

for i in {1..22}; do

sbatch -p highmem -N 1 -n 16 heightA.sh $i;

done
