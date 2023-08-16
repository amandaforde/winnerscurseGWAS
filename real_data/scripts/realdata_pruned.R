## WINNER'S CURSE REAL DATA SCRIPT 2 - PRUNED DATA SETS:

## This script is very similar to the first script, 'winnerscurse_realdata.R'.
## However, in this script, pruned versions of the six real data sets are first
## obtained. Following this, the winner's curse correction methods are then
## evaluated at two significance thresholds, 5e-8 and 5e-4, for each pruned data
## set, using the estimated MSE among significant SNPs.

## Load required packages:
library(winnerscurse)
library(ggplot2)
library(dplyr)
library(scam)
library(tidyr)
library("RColorBrewer")
col <- brewer.pal(8,"Dark2")
col1 <- brewer.pal(11,"RdYlBu")

set.seed(1998)

################################################################################
################################################################################

## load set of pruned SNPs and obtain pruned data sets
pruned_SNPs <- read.table('real_data/data/pruned_SNPs_bmi_1.txt',header=TRUE)

summary_data_bmi_1 <-  read.table('real_data/data/summary_data_bmi_1.txt',header=TRUE)
summary_data_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$rsid %in% pruned_SNPs$V1,]

summary_data_bmi_2 <-  read.table('real_data/data/summary_data_bmi_2.txt',header=TRUE)
summary_data_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$rsid %in% pruned_SNPs$V1,]

summary_data_T2D_1 <-  read.table('real_data/data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)
summary_data_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$rsid %in% pruned_SNPs$V1,]

summary_data_T2D_2 <-  read.table('real_data/data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)
summary_data_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$rsid %in% pruned_SNPs$V1,]

summary_data_height_1 <-  read.table('real_data/data/summary_data_height_1.txt',header=TRUE)
summary_data_height_1 <- summary_data_height_1[summary_data_height_1$rsid %in% pruned_SNPs$V1,]

summary_data_height_2 <-  read.table('real_data/data/summary_data_height_2.txt',header=TRUE)
summary_data_height_2 <- summary_data_height_2[summary_data_height_2$rsid %in% pruned_SNPs$V1,]

################################################################################

## 1) BMI 1:
## COMPUTING ESTIMATED MSE FOR PRUNED BMI 1 AT THRESHOLDS 5e-8 AND 5e-4

# Apply methods
summary_stats <- summary_data_bmi_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

## MSE (5e-8)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_bmi_1$beta_rep

mse_BMI1 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## MSE (5e-4)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4) ## need to reapply conditional likelihood
summary_data_sig2 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

mse_BMI1_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

################################################################################

## 2) BMI 2:
## COMPUTING ESTIMATED MSE FOR PRUNED BMI 2 AT THRESHOLDS 5e-8 AND 5e-4

# Apply methods
summary_stats <- summary_data_bmi_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

## MSE (5e-8)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_bmi_2$beta_rep

mse_BMI2 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## MSE (5e-4)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4) ## need to reapply conditional likelihood
summary_data_sig2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

mse_BMI2_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

################################################################################

## 3) T2D 1:
## COMPUTING ESTIMATED MSE FOR PRUNED T2D 1 AT THRESHOLDS 5e-8 AND 5e-4

# Apply methods
summary_stats <- summary_data_T2D_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

## MSE (5e-8)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]
summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_T2D_1$beta_rep

mse_T2D1 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## MSE (5e-4)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4) ## need to reapply conditional likelihood
summary_data_sig2 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

mse_T2D1_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

################################################################################

## 4) T2D 2:
## COMPUTING ESTIMATED MSE AND BIAS FOR T2D 2 AT THRESHOLDS 5e-8 AND 5e-4

# Apply methods
summary_stats <- summary_data_T2D_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

## MSE (5e-8)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]
summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_T2D_2$beta_rep

mse_T2D2 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## MSE (5e-4)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4) ## need to reapply conditional likelihood
summary_data_sig2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

mse_T2D2_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

###############################################################################

## 5) Height 1:
## COMPUTING ESTIMATED MSE AND BIAS FOR Height 1 AT THRESHOLDS 5e-8 AND 5e-4

# Apply methods
summary_stats <- summary_data_height_1[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

## MSE (5e-8)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_height_1$beta_rep

mse_height1 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## MSE (5e-4)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4) ## need to reapply conditional likelihood
summary_data_sig2 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

mse_height1_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

################################################################################

## 6) Height 2:
## COMPUTING ESTIMATED MSE AND BIAS FOR Height 2 AT THRESHOLDS 5e-8 AND 5e-4

# Apply methods
summary_stats <- summary_data_height_2[,3:5]
out_CL <- conditional_likelihood(summary_data = summary_stats)
out_EB <- empirical_bayes(summary_data = summary_stats)
out_EB_df <- empirical_bayes(summary_data=summary_stats, method="fix_df")
out_EB_scam <- empirical_bayes(summary_data=summary_stats,method="scam")
out_EB_gam_po <- empirical_bayes(summary_data=summary_stats,method="gam_po")
out_EB_gam_nb <- empirical_bayes(summary_data=summary_stats,method="gam_nb")
out_FIQT <- FDR_IQT(summary_data = summary_stats)
out_BR <- BR_ss(summary_data = summary_stats)

## MSE (5e-8)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-8,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-8,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-8,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-8,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-8,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-8,]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-8,]
true_beta <- summary_data_sig_height_2$beta_rep

mse_height2 <- data.frame(naive = mean((true_beta - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

## MSE (5e-4)
out_EB_sig <- out_EB[2*(pnorm(abs(out_EB$beta/out_EB$se), lower.tail=FALSE)) < 5e-4,]
out_EB_df_sig <- out_EB_df[2*(pnorm(abs(out_EB_df$beta/out_EB_df$se), lower.tail=FALSE)) < 5e-4,]
out_EB_scam_sig <- out_EB_scam[2*(pnorm(abs(out_EB_scam$beta/out_EB_scam$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gampo_sig <- out_EB_gam_po[2*(pnorm(abs(out_EB_gam_po$beta/out_EB_gam_po$se), lower.tail=FALSE)) < 5e-4,]
out_EB_gamnb_sig <- out_EB_gam_nb[2*(pnorm(abs(out_EB_gam_nb$beta/out_EB_gam_nb$se), lower.tail=FALSE)) < 5e-4,]
out_FIQT_sig <- out_FIQT[2*(pnorm(abs(out_FIQT$beta/out_FIQT$se), lower.tail=FALSE)) < 5e-4,]
out_BR_ss_sig <- out_BR[2*(pnorm(abs(out_BR$beta/out_BR$se), lower.tail=FALSE)) < 5e-4,]
out_CL <- conditional_likelihood(summary_data = summary_stats, alpha=5e-4) ## need to reapply conditional likelihood
summary_data_sig2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-4,]
true_beta2 <- summary_data_sig2$beta_rep

mse_height2_5e_4 <- data.frame(naive = mean((true_beta2 - out_EB_sig$beta)^2) - mean(out_EB_sig$se^2), CL1 = mean((true_beta2 - out_CL$beta.cl1)^2) - mean(out_EB_sig$se^2), CL2 = mean((true_beta2 - out_CL$beta.cl2)^2) - mean(out_EB_sig$se^2), CL3 = mean((true_beta2 - out_CL$beta.cl3)^2) - mean(out_EB_sig$se^2), EB = mean((true_beta2 - out_EB_sig$beta_EB)^2) - mean(out_EB_sig$se^2),  EB_df = mean((true_beta2 - out_EB_df_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_scam = mean((true_beta2 - out_EB_scam_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_po = mean((true_beta2 - out_EB_gampo_sig$beta_EB)^2) - mean(out_EB_sig$se^2), EB_gam_nb = mean((true_beta2 - out_EB_gamnb_sig$beta_EB)^2) - mean(out_EB_sig$se^2), boot = mean((true_beta2 - out_BR_ss_sig$beta_BR_ss)^2) - mean(out_EB_sig$se^2), FIQT = mean((true_beta2 - out_FIQT_sig$beta_FIQT)^2) - mean(out_EB_sig$se^2))

################################################################################

## Combine above results and save them in results file

# S15 Table: PRUNED MSE 5e-8
mse_5e_8 <- rbind(mse_BMI1,mse_BMI2, mse_T2D1, mse_T2D2, mse_height1, mse_height2)
mse_5e_8 <- round(mse_5e_8[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_8 <- cbind(GWAS,mse_5e_8)
write.csv(mse_5e_8,"real_data/results/pruned_mse_5e_8.txt", row.names = FALSE)

# S16 Table: PRUNED MSE 5e-4
mse_5e_4 <- rbind(mse_BMI1_5e_4,mse_BMI2_5e_4, mse_T2D1_5e_4, mse_T2D2_5e_4, mse_height1_5e_4, mse_height2_5e_4)
mse_5e_4 <- round(mse_5e_4[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_4 <- cbind(GWAS,mse_5e_4)
write.csv(mse_5e_4,"real_data/results/pruned_mse_5e_4.txt", row.names = FALSE)

################################################################################

