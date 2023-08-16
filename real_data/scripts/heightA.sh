#!/bin/sh

input=$1

plink2 --pfile /data/lfahey/ukb_qc/qcd_files2/${input}_qcd --rm-dup force-first --keep height_gwasA.txt --pheno height_gwasA.txt --make-bed --out /data/aforde/height_gwasA/${input}_qcd_heightA

plink2 --bfile /data/aforde/height_gwasA/${input}_qcd_heightA --glm omit-ref hide-covar --covar covariates.txt --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 AGE --out /data/aforde/results_height/heightA_res_${input}
