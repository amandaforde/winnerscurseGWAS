#!/bin/sh

input=$1


plink2 --pfile /data/lfahey/ukb_qc/qcd_files2/${input}_qcd --rm-dup force-first --keep bmi_gwasB.txt --pheno bmi_gwasB.txt --make-bed --out /data/aforde/bmi_gwasB/${input}_qcd_bmiB

plink2 --bfile /data/aforde/bmi_gwasB/${input}_qcd_bmiB --glm omit-red hide-covar --covar covariates.txt --covar-variance-standardize PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 AGE --out /data/aforde/results_bmi/bmiB_res_${input}
