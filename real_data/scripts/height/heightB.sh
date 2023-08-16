#!/bin/sh

input=$1


plink2 --pfile /data/lfahey/ukb_qc/qcd_files2/${input}_qcd --rm-dup force-first --keep height_gwasB.txt --pheno height_gwasB.txt --make-bed --out /data/aforde/height_gwasB/${input}_qcd_heightB

plink2 --bfile /data/aforde/height_gwasB/${input}_qcd_heightB --glm omit-red hide-covar --covar covariates.txt --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 AGE --out /data/aforde/results_height/heightB_res_${input}
