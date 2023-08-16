#!/bin/sh -l
#SBATCH -J bmi_2
#SBATCH --partition=highmem
#SBATCH -N 1
#SBATCH -n 1

Rscript summary_data_bmi.R results_bmi/bmiB_res_1.PHENO1.glm.linear results_bmi/bmiA_res_1.PHENO1.glm.linear results_bmi/bmiB_res_2.PHENO1.glm.linear results_bmi/bmiA_res_2.PHENO1.glm.linear results_bmi/bmiB_res_3.PHENO1.glm.linear results_bmi/bmiA_res_3.PHENO1.glm.linear results_bmi/bmiB_res_4.PHENO1.glm.linear results_bmi/bmiA_res_4.PHENO1.glm.linear results_bmi/bmiB_res_5.PHENO1.glm.linear results_bmi/bmiA_res_5.PHENO1.glm.linear results_bmi/bmiB_res_6.PHENO1.glm.linear results_bmi/bmiA_res_6.PHENO1.glm.linear results_bmi/bmiB_res_7.PHENO1.glm.linear results_bmi/bmiA_res_7.PHENO1.glm.linear results_bmi/bmiB_res_8.PHENO1.glm.linear results_bmi/bmiA_res_8.PHENO1.glm.linear results_bmi/bmiB_res_9.PHENO1.glm.linear results_bmi/bmiA_res_9.PHENO1.glm.linear results_bmi/bmiB_res_10.PHENO1.glm.linear results_bmi/bmiA_res_10.PHENO1.glm.linear results_bmi/bmiB_res_11.PHENO1.glm.linear results_bmi/bmiA_res_11.PHENO1.glm.linear results_bmi/bmiB_res_12.PHENO1.glm.linear results_bmi/bmiA_res_12.PHENO1.glm.linear results_bmi/bmiB_res_13.PHENO1.glm.linear results_bmi/bmiA_res_13.PHENO1.glm.linear results_bmi/bmiB_res_14.PHENO1.glm.linear results_bmi/bmiA_res_14.PHENO1.glm.linear results_bmi/bmiB_res_15.PHENO1.glm.linear results_bmi/bmiA_res_15.PHENO1.glm.linear results_bmi/bmiB_res_16.PHENO1.glm.linear results_bmi/bmiA_res_16.PHENO1.glm.linear results_bmi/bmiB_res_17.PHENO1.glm.linear results_bmi/bmiA_res_17.PHENO1.glm.linear results_bmi/bmiB_res_18.PHENO1.glm.linear results_bmi/bmiA_res_18.PHENO1.glm.linear results_bmi/bmiB_res_19.PHENO1.glm.linear results_bmi/bmiA_res_19.PHENO1.glm.linear results_bmi/bmiB_res_20.PHENO1.glm.linear results_bmi/bmiA_res_20.PHENO1.glm.linear results_bmi/bmiB_res_21.PHENO1.glm.linear results_bmi/bmiA_res_21.PHENO1.glm.linear results_bmi/bmiB_res_22.PHENO1.glm.linear results_bmi/bmiA_res_22.PHENO1.glm.linear summary_data_bmi_2.txt
