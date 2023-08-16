#!/bin/sh
#SBATCH -J ld_bmi
#SBATCH --partition=highmem
#SBATCH -N 1
#SBATCH -n 16

for i in {1..22}; do plink2 --bfile /data/aforde/bmi_gwasA/${i}_qcd_bmiA --indep-pairwise 50 5 0.5 --out ${i}_bmiA_ld; done
 
Rscript join_pruned.R 1_bmiA_ld.prune.in 2_bmiA_ld.prune.in 3_bmiA_ld.prune.in 4_bmiA_ld.prune.in 5_bmiA_ld.prune.in 6_bmiA_ld.prune.in 7_bmiA_ld.prune.in 8_bmiA_ld.prune.in 9_bmiA_ld.prune.in 10_bmiA_ld.prune.in 11_bmiA_ld.prune.in 12_bmiA_ld.prune.in 13_bmiA_ld.prune.in 14_bmiA_ld.prune.in 15_bmiA_ld.prune.in 16_bmiA_ld.prune.in 17_bmiA_ld.prune.in 18_bmiA_ld.prune.in 19_bmiA_ld.prune.in 20_bmiA_ld.prune.in 21_bmiA_ld.prune.in 22_bmiA_ld.prune.in pruned_SNPS_bmi_1.txt 
