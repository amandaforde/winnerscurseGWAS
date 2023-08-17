## WINNER'S CURSE REAL DATA SCRIPT 3 - INDEPENDENT SIGNALS:

## This script accompanies our attempt at determining the number of independent
## signals required to ensure appropriate performance of the winner's curse
## correction methods. Firstly, the estimated MSE is computed for the BMI data
## sets at thresholds 5e-10, 5e-12 and 5e-14, and for the height data sets at
## thresholds 5e-32, 5e-34 and 5e-36. These results are illustrated. Manhattan
## plots for each data set are also produced. 

## Outputs:
## 1. S23_Fig.tiff
## 2. S24_Fig.tiff
## 3. S25_Fig.tiff
## 4. S26_Fig.tiff
## 5. S27_Fig.tiff
## 6. S28_Fig.tiff 

## Load required packages:
library(winnerscurse)
library(ggplot2)
library(dplyr)
library(scam)
library(tidyr)
library(patchwork)
library(ggpubr)
library(qqman)
library("RColorBrewer")
col <- brewer.pal(8,"Dark2")
col1 <- brewer.pal(11,"RdYlBu")

set.seed(1998)

################################################################################
################################################################################
## PART A)

## 1) BMI: 
## COMPUTING ESTIMATED MSE AND BIAS FOR BMI 1 AND 2 AT THRESHOLDS 5e-10,5e-12
## AND 5e-14

summary_data_bmi_1 <- read.table('real_data/data/summary_data_bmi_1.txt',header=TRUE)
summary_data_bmi_2 <- read.table('real_data/data/summary_data_bmi_2.txt',header=TRUE)

## BMI 1
summary_stats <- summary_data_bmi_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-10)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# MSE 5e-10
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-10,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-10,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-10,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-10,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-10,]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-10,]
true_beta <- summary_data_sig_bmi_1$beta_rep
mse_bmi1_5e_10 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-12
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
mse_bmi1_5e_12 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-14
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-14,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-14,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-14,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-14,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-14,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-14)
summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-14,]
true_beta2 <- summary_data_sig2$beta_rep
mse_bmi1_5e_14 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## BMI 2

summary_stats <- summary_data_bmi_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-10)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# MSE 5e-10
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-10,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-10,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-10,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-10,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-10,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-10,]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-10,]
true_beta <- summary_data_sig_bmi_2$beta_rep
mse_bmi2_5e_10 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-12
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
mse_bmi2_5e_12 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-14
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-14,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-14,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-14,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-14,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-14,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-14,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-14)
summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-14,]
true_beta2 <- summary_data_sig2$beta_rep
mse_bmi2_5e_14 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Combine all results 

mse_5e_10 <- rbind(mse_bmi1_5e_10, mse_bmi2_5e_10)
mse_5e_10 <- round(mse_5e_10[,1:11],5)
GWAS <- c("BMI 1", "BMI 2")
mse_5e_10 <- cbind(GWAS,mse_5e_10)

mse_5e_12 <- rbind(mse_bmi1_5e_12, mse_bmi2_5e_12)
mse_5e_12 <- round(mse_5e_12[,1:11],5)
GWAS <- c("BMI 1", "BMI 2")
mse_5e_12 <- cbind(GWAS,mse_5e_12)

mse_5e_14 <- rbind(mse_bmi1_5e_14, mse_bmi2_5e_14)
mse_5e_14 <- round(mse_5e_14[,1:11],5)
GWAS <- c("BMI 1", "BMI 2")
mse_5e_14 <- cbind(GWAS,mse_5e_14)

## Plot results 

new_mse_5e_10 <- pivot_longer(mse_5e_10, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_12 <- pivot_longer(mse_5e_12, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_14 <- pivot_longer(mse_5e_14, -c(GWAS), values_to = "MSE", names_to = "Method")

new_mse_5e_10$Method <- factor(new_mse_5e_10$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_10$GWAS <- factor(new_mse_5e_10$GWAS, levels=c("BMI 1", "BMI 2"))
new_mse_5e_12$Method <- factor(new_mse_5e_12$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_12$GWAS <- factor(new_mse_5e_12$GWAS, levels=c("BMI 1", "BMI 2"))
new_mse_5e_14$Method <- factor(new_mse_5e_14$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_14$GWAS <- factor(new_mse_5e_14$GWAS, levels=c("BMI 1", "BMI 2"))

nsig_BMI1_5e_10 <- sum(2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-10)
nsig_BMI1_5e_12 <-sum(2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-12)
nsig_BMI1_5e_14 <-sum(2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-14)
nsig_BMI2_5e_10 <- sum(2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-10)
nsig_BMI2_5e_12 <- sum(2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-12)
nsig_BMI2_5e_14 <-sum(2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-14)

ann_text1 <- data.frame(MSE = 0.00145,Method = "boot",lab = "label", GWAS = factor("BMI 1",levels = c("BMI 1", "BMI 2")))
ann_text2 <- data.frame(MSE = 0.00145,Method = "boot",lab = "label", GWAS = factor("BMI 2",levels = c("BMI 1", "BMI 2")))
plot_10 <- ggplot(new_mse_5e_10,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-10))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_10[new_mse_5e_10$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) +
  geom_label(data = ann_text1,label = "3333 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "4107 sig SNPs",size=3.5,fill="#F4F4F4")

ann_text1 <- data.frame(MSE = 0.0019,Method = "boot",lab = "label", GWAS = factor("BMI 1",levels = c("BMI 1", "BMI 2")))
ann_text2 <- data.frame(MSE = 0.0019,Method = "boot",lab = "label", GWAS = factor("BMI 2",levels = c("BMI 1", "BMI 2")))
plot_12 <- ggplot(new_mse_5e_12,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-12))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_12[new_mse_5e_12$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) +
  geom_label(data = ann_text1,label = "2113 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "2745 sig SNPs",size=3.5,fill="#F4F4F4")

ann_text1 <- data.frame(MSE = 0.002,Method = "boot",lab = "label", GWAS = factor("BMI 1",levels = c("BMI 1", "BMI 2")))
ann_text2 <- data.frame(MSE = 0.002,Method = "boot",lab = "label", GWAS = factor("BMI 2",levels = c("BMI 1", "BMI 2")))
plot_14 <- ggplot(new_mse_5e_14,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-14))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_14[new_mse_5e_14$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) +
  geom_label(data = ann_text1,label = "1519 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "1742 sig SNPs",size=3.5,fill="#F4F4F4")

## S23 Fig. Estimated MSE of significant SNPs at thresholds 5 × 10-10, 5 × 10-12
## and 5 × 10-14 for each method and BMI data set.
## Save as 'figures/S23_Fig.tiff' with dimensions 900 x 1200 pixels.

figure <- plot_10 + plot_12 + plot_14
figure + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

##################################################################################
##################################################################################

## 2) Height: 
## COMPUTING ESTIMATED MSE AND BIAS FOR Height 1 AND 2 AT THRESHOLDS 5e-32,5e-34
## AND 5e-36

summary_data_height_1 <- read.table('real_data/data/summary_data_height_1.txt',header=TRUE)
summary_data_height_2 <- read.table('real_data/data/summary_data_height_2.txt',header=TRUE)

## Height 1
summary_stats <- summary_data_height_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-32)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# MSE 5e-32
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-32,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-32,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-32,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-32,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-32,]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-32,]
true_beta <- summary_data_sig_height_1$beta_rep
mse_height1_5e_32 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-34
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-34,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-34,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-34,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-34,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-34,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-34,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-34,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-34)
summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-34,]
true_beta2 <- summary_data_sig2$beta_rep
mse_height1_5e_34 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-36
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-36,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-36,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-36,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-36,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-36,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-36)
summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-36,]
true_beta2 <- summary_data_sig2$beta_rep
mse_height1_5e_36 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## Height 2
summary_stats <- summary_data_height_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-32)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

# MSE 5e-32
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-32,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-32,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-32,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-32,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-32,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-32,]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-32,]
true_beta <- summary_data_sig_height_2$beta_rep
mse_height2_5e_32 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-34
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-12,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-12,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-12,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-12,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-12,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-12,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-12)
summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-12,]
true_beta2 <- summary_data_sig2$beta_rep
mse_height2_5e_34 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

# MSE 5e-36
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-36,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-36,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-36,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-36,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-36,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-36,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-36)
summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-36,]
true_beta2 <- summary_data_sig2$beta_rep
mse_height2_5e_36 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## Combine all results
mse_5e_32 <- rbind(mse_height1_5e_32, mse_height2_5e_32)
mse_5e_32 <- round(mse_5e_32[,1:11],5)
GWAS <- c("Height 1", "Height 2")
mse_5e_32 <- cbind(GWAS,mse_5e_32)

mse_5e_34 <- rbind(mse_height1_5e_34, mse_height2_5e_34)
mse_5e_34 <- round(mse_5e_34[,1:11],5)
GWAS <- c("Height 1", "Height 2")
mse_5e_34 <- cbind(GWAS,mse_5e_34)

mse_5e_36 <- rbind(mse_height1_5e_36, mse_height2_5e_36)
mse_5e_36 <- round(mse_5e_36[,1:11],5)
GWAS <- c("Height 1", "Height 2")
mse_5e_36 <- cbind(GWAS,mse_5e_36)

## Plot results

new_mse_5e_32 <- pivot_longer(mse_5e_32, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_34 <- pivot_longer(mse_5e_34, -c(GWAS), values_to = "MSE", names_to = "Method")
new_mse_5e_36 <- pivot_longer(mse_5e_36, -c(GWAS), values_to = "MSE", names_to = "Method")

new_mse_5e_32$Method <- factor(new_mse_5e_32$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_32$GWAS <- factor(new_mse_5e_32$GWAS, levels=c("Height 1", "Height 2"))
new_mse_5e_34$Method <- factor(new_mse_5e_34$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_34$GWAS <- factor(new_mse_5e_34$GWAS, levels=c("Height 1", "Height 2"))
new_mse_5e_36$Method <- factor(new_mse_5e_36$Method, levels=c("naive", "CL1", "CL2", "CL3", "EB", "EB_df", "EB_gam_nb", "EB_gam_po", "EB_scam", "boot", "FIQT"))
new_mse_5e_36$GWAS <- factor(new_mse_5e_36$GWAS, levels=c("Height 1", "Height 2"))

nsig_height1_5e_32 <- sum(2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-32)
nsig_height1_5e_34 <- sum(2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-34)
nsig_height1_5e_36 <- sum(2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-36)
nsig_height2_5e_32 <- sum(2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-32)
nsig_height2_5e_34 <- sum(2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-34)
nsig_height2_5e_36 <- sum(2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-36)

ann_text1 <- data.frame(MSE = 0.0062,Method = "boot",lab = "label", GWAS = factor("Height 1",levels = c("Height 1", "Height 2")))
ann_text2 <- data.frame(MSE = 0.0062,Method = "boot",lab = "label", GWAS = factor("Height 2",levels = c("Height 1", "Height 2")))
plot_32 <- ggplot(new_mse_5e_32,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(.~GWAS,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-32))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle = 0)) + geom_vline(data = new_mse_5e_32[new_mse_5e_32$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) +
  geom_label(data = ann_text1,label = "3459 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "4163 sig SNPs",size=3.5,fill="#F4F4F4")

ann_text1 <- data.frame(MSE = 0.0067,Method = "boot",lab = "label", GWAS = factor("Height 1",levels = c("Height 1", "Height 2")))
ann_text2 <- data.frame(MSE = 0.0067,Method = "boot",lab = "label", GWAS = factor("Height 2",levels = c("Height 1", "Height 2")))
plot_34 <- ggplot(new_mse_5e_34,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-34))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_34[new_mse_5e_34$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) +
  geom_label(data = ann_text1,label = "3062 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "3358 sig SNPs",size=3.5,fill="#F4F4F4")

ann_text1 <- data.frame(MSE = 0.0047,Method = "boot",lab = "label", GWAS = factor("Height 1",levels = c("Height 1", "Height 2")))
ann_text2 <- data.frame(MSE = 0.0047,Method = "boot",lab = "label", GWAS = factor("Height 2",levels = c("Height 1", "Height 2")))
plot_36 <- ggplot(new_mse_5e_36,aes(x=MSE,y=Method,fill=Method)) + geom_bar(stat="identity")+ geom_col(width=1) + facet_wrap(GWAS~.,nrow=6, strip.position = "right") + scale_shape_manual(values=c(6,20,20,20,15,18,18,18,18,10,15)) + scale_fill_manual(values=c(col[1],col[2],col[8],col[4],col[5],col[3],col[6],col[7],col1[1],col1[11],col1[9])) + ylab("Method") +
  xlab(expression(paste(italic("Estimated"), " MSE at ", 5%*%10^-36))) + theme_bw() + theme(text = element_text(size=12),legend.position = "none",panel.grid.major = element_blank(), panel.grid.minor = element_blank(), strip.text.y = element_text(angle=0)) + geom_vline(data = new_mse_5e_36[new_mse_5e_36$Method=="naive",], aes(xintercept = MSE), linetype=2, size=0.5) + geom_vline(xintercept=0) +
  geom_label(data = ann_text1,label = "2748 sig SNPs",size=3.5,fill="#F4F4F4") +
  geom_label(data = ann_text2,label = "2867 sig SNPs",size=3.5,fill="#F4F4F4")

## S24 Fig. Estimated MSE of significant SNPs at thresholds 5 × 10-10, 5 × 10-12
## and 5 × 10-14 for each method and BMI data set.
## Save as 'figures/S24_Fig.tiff' with dimensions 900 x 1200 pixels.

figure <- plot_32 + plot_34 + plot_36
figure + plot_layout(guides = "collect", ncol=1) + plot_annotation(tag_levels = 'A') &
  theme(plot.tag = element_text(face = "bold"))

################################################################################
################################################################################

## PART B) Creating Manhattan plots for BMI 1, BMI 2, Height 1 and Height 2

summary_data_bmi_1_sub <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 0.05,]
summary_data_bmi_1_new <- data.frame(SNP = summary_data_bmi_1_sub$rsid, CHR = summary_data_bmi_1_sub$chr, BP = summary_data_bmi_1_sub$pos, P = 2*(pnorm(abs(summary_data_bmi_1_sub$beta/summary_data_bmi_1_sub$se), lower.tail=FALSE)))

## S25 Fig. Manhattan plot for the first BMI data set. 
## Save as 'figures/S25_Fig.tiff' with dimensions 1200 x 600 pixels.
manhattan(summary_data_bmi_1_new, main = "BMI 1", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-10))


summary_data_bmi_2_sub <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 0.05,]
summary_data_bmi_2_new <- data.frame(SNP = summary_data_bmi_2_sub$rsid, CHR = summary_data_bmi_2_sub$chr, BP = summary_data_bmi_2_sub$pos, P = 2*(pnorm(abs(summary_data_bmi_2_sub$beta/summary_data_bmi_2_sub$se), lower.tail=FALSE)))

## S26 Fig. Manhattan plot for the second BMI data set. 
## Save as 'figures/S26_Fig.tiff' with dimensions 1200 x 600 pixels.
manhattan(summary_data_bmi_2_new, main = "BMI 2", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-12))


summary_data_height_1_sub <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 0.05,]
summary_data_height_1_new <- data.frame(SNP = summary_data_height_1_sub$rsid, CHR = summary_data_height_1_sub$chr, BP = summary_data_height_1_sub$pos, P = 2*(pnorm(abs(summary_data_height_1_sub$beta/summary_data_height_1_sub$se), lower.tail=FALSE)))

## S27 Fig. Manhattan plot for the first height data set. 
## Save as 'figures/S27_Fig.tiff' with dimensions 1200 x 600 pixels.
manhattan(summary_data_height_1_new, main = "Height 1", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-32))


summary_data_height_2_sub <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 0.05,]
summary_data_height_2_new <- data.frame(SNP = summary_data_height_2_sub$rsid, CHR = summary_data_height_2_sub$chr, BP = summary_data_height_2_sub$pos, P = 2*(pnorm(abs(summary_data_height_2_sub$beta/summary_data_height_2_sub$se), lower.tail=FALSE)))

## S28Fig. Manhattan plot for the second height data set. 
## Save as 'figures/S28_Fig.tiff' with dimensions 1200 x 600 pixels.
manhattan(summary_data_height_2_new, main = "Height 2", cex = 0.6, cex.axis = 0.9, suggestiveline = -log10(5e-34))

