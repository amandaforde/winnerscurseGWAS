#!/bin/sh

input=$1

plink2 --pfile /data/lfahey/ukb_qc/qcd_files2/${input}_qcd --rm-dup force-first --keep T2D_gwasA.txt --pheno T2D_gwasA.txt --make-bed --out /data/aforde/T2D_gwasA/${input}_qcd_T2DA

plink2 --bfile /data/aforde/T2D_gwasA/${input}_qcd_T2DA --glm firth-fallback omit-ref hide-covar --covar covariates.txt --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 AGE --out /data/aforde/results_T2D/T2DA_res_${input}


