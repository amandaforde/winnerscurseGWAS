## WINNER'S CURSE REAL DATA SCRIPT 1:

## This script uses 6 sets of GWAS summary statistics, two related to each
## trait. The traits of interest are body mass index (T2D), type 2 diabetes
## (T2D) and height. Each data set contains six columns, named 'chr', 'pos',
## 'rsid', 'beta', 'se' and 'beta_rep'. Here 'beta_rep' contains the association
## estimate obtained for each SNP in the corresponding replication GWAS.

## This script's primary purpose is to evaluate a number of winner's curse
## correction methods by comparing their performance at two significance
## thresholds, 5e-8 and 5e-4, for each data set, using both the estimated MSE
## among significant SNPs and the average bias over all significant SNPs.

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
## PART A)

## In the first part, we perform an initial exploration of the data sets. The
## number of significant SNPs at the thresholds, 5e-8 and 5e-4, as well as
## proportions that indicate the extent of winner's curse for each data set are
## computed. Firstly, we determine the proportion of these significant SNPs that
## had smaller estimated effect sizes, in terms of absolute value, in their
## respective replication GWAS, and secondly, the proportion of significant SNPs
## that have significantly overestimated effect sizes.

## BMI 1
summary_data_bmi_1 <-  read.table('real_data/data/summary_data_bmi_1.txt',header=TRUE)

sig_BMI1 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-8,]
nsig_BMI1 <- nrow(sig_BMI1)
prop_over_BMI1 <- sum(abs(sig_BMI1$beta)>abs(sig_BMI1$beta_rep))/nrow(sig_BMI1)
prop_sigover_BMI1 <- sum(abs(sig_BMI1$beta)>(abs(sig_BMI1$beta_rep)+1.96*sig_BMI1$se))/nrow(sig_BMI1)

sig_BMI1_5e_4 <- summary_data_bmi_1[2*(pnorm(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se), lower.tail=FALSE)) < 5e-4,]
nsig_BMI1_5e_4 <- nrow(sig_BMI1_5e_4)
prop_over_BMI1_5e_4 <- sum(abs(sig_BMI1_5e_4$beta)>abs(sig_BMI1_5e_4$beta_rep))/nrow(sig_BMI1_5e_4)
prop_sigover_BMI1_5e_4 <- sum(abs(sig_BMI1_5e_4$beta)>(abs(sig_BMI1_5e_4$beta_rep)+1.96*sig_BMI1_5e_4$se))/nrow(sig_BMI1_5e_4)

values_BMI1 <- c(nsig_BMI1,prop_over_BMI1,prop_sigover_BMI1,nsig_BMI1_5e_4,prop_over_BMI1_5e_4,prop_sigover_BMI1_5e_4)

## BMI 2
summary_data_bmi_2 <-  read.table('real_data/data/summary_data_bmi_2.txt',header=TRUE)

sig_BMI2 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-8,]
nsig_BMI2 <- nrow(sig_BMI2)
prop_over_BMI2 <- sum(abs(sig_BMI2$beta)>abs(sig_BMI2$beta_rep))/nrow(sig_BMI2)
prop_sigover_BMI2 <- sum(abs(sig_BMI2$beta)>(abs(sig_BMI2$beta_rep)+1.96*sig_BMI2$se))/nrow(sig_BMI2)

sig_BMI2_5e_4 <- summary_data_bmi_2[2*(pnorm(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se), lower.tail=FALSE)) < 5e-4,]
nsig_BMI2_5e_4 <- nrow(sig_BMI2_5e_4)
prop_over_BMI2_5e_4 <- sum(abs(sig_BMI2_5e_4$beta)>abs(sig_BMI2_5e_4$beta_rep))/nrow(sig_BMI2_5e_4)
prop_sigover_BMI2_5e_4 <- sum(abs(sig_BMI2_5e_4$beta)>(abs(sig_BMI2_5e_4$beta_rep)+1.96*sig_BMI2_5e_4$se))/nrow(sig_BMI2_5e_4)

values_BMI2 <- c(nsig_BMI2,prop_over_BMI2,prop_sigover_BMI2,nsig_BMI2_5e_4,prop_over_BMI2_5e_4,prop_sigover_BMI2_5e_4)

## T2D 1
summary_data_T2D_1 <-  read.table('real_data/data/summary_data_T2D_1.txt',header=TRUE)
summary_data_T2D_1$beta <- log(summary_data_T2D_1$beta)
summary_data_T2D_1$beta_rep <- log(summary_data_T2D_1$beta_rep)

sig_T2D1 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-8,]
nsig_T2D1 <- nrow(sig_T2D1)
prop_over_T2D1 <- sum(abs(sig_T2D1$beta)>abs(sig_T2D1$beta_rep))/nrow(sig_T2D1)
prop_sigover_T2D1 <- sum(abs(sig_T2D1$beta)>(abs(sig_T2D1$beta_rep)+1.96*sig_T2D1$se))/nrow(sig_T2D1)

sig_T2D1_5e_4 <- summary_data_T2D_1[2*(pnorm(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se), lower.tail=FALSE)) < 5e-4,]
nsig_T2D1_5e_4 <- nrow(sig_T2D1_5e_4)
prop_over_T2D1_5e_4 <- sum(abs(sig_T2D1_5e_4$beta)>abs(sig_T2D1_5e_4$beta_rep))/nrow(sig_T2D1_5e_4)
prop_sigover_T2D1_5e_4 <- sum(abs(sig_T2D1_5e_4$beta)>(abs(sig_T2D1_5e_4$beta_rep)+1.96*sig_T2D1_5e_4$se))/nrow(sig_T2D1_5e_4)

values_T2D1 <- c(nsig_T2D1,prop_over_T2D1,prop_sigover_T2D1,nsig_T2D1_5e_4,prop_over_T2D1_5e_4,prop_sigover_T2D1_5e_4)

## T2D 2
summary_data_T2D_2 <-  read.table('real_data/data/summary_data_T2D_2.txt',header=TRUE)
summary_data_T2D_2$beta <- log(summary_data_T2D_2$beta)
summary_data_T2D_2$beta_rep <- log(summary_data_T2D_2$beta_rep)

sig_T2D2 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-8,]
nsig_T2D2 <- nrow(sig_T2D2)
prop_over_T2D2 <- sum(abs(sig_T2D2$beta)>abs(sig_T2D2$beta_rep))/nrow(sig_T2D2)
prop_sigover_T2D2 <- sum(abs(sig_T2D2$beta)>(abs(sig_T2D2$beta_rep)+1.96*sig_T2D2$se))/nrow(sig_T2D2)

sig_T2D2_5e_4 <- summary_data_T2D_2[2*(pnorm(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se), lower.tail=FALSE)) < 5e-4,]
nsig_T2D2_5e_4 <- nrow(sig_T2D2_5e_4)
prop_over_T2D2_5e_4 <- sum(abs(sig_T2D2_5e_4$beta)>abs(sig_T2D2_5e_4$beta_rep))/nrow(sig_T2D2_5e_4)
prop_sigover_T2D2_5e_4 <- sum(abs(sig_T2D2_5e_4$beta)>(abs(sig_T2D2_5e_4$beta_rep)+1.96*sig_T2D2_5e_4$se))/nrow(sig_T2D2_5e_4)

values_T2D2 <- c(nsig_T2D2,prop_over_T2D2,prop_sigover_T2D2,nsig_T2D2_5e_4,prop_over_T2D2_5e_4,prop_sigover_T2D2_5e_4)

## Height 1
summary_data_height_1 <-  read.table('real_data/data/summary_data_height_1.txt',header=TRUE)

sig_height1 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-8,]
nsig_height1 <- nrow(sig_height1)
prop_over_height1 <- sum(abs(sig_height1$beta)>abs(sig_height1$beta_rep))/nrow(sig_height1)
prop_sigover_height1 <- sum(abs(sig_height1$beta)>(abs(sig_height1$beta_rep)+1.96*sig_height1$se))/nrow(sig_height1)

sig_height1_5e_4 <- summary_data_height_1[2*(pnorm(abs(summary_data_height_1$beta/summary_data_height_1$se), lower.tail=FALSE)) < 5e-4,]
nsig_height1_5e_4 <- nrow(sig_height1_5e_4)
prop_over_height1_5e_4 <- sum(abs(sig_height1_5e_4$beta)>abs(sig_height1_5e_4$beta_rep))/nrow(sig_height1_5e_4)
prop_sigover_height1_5e_4 <- sum(abs(sig_height1_5e_4$beta)>(abs(sig_height1_5e_4$beta_rep)+1.96*sig_height1_5e_4$se))/nrow(sig_height1_5e_4)

values_height1 <- c(nsig_height1,prop_over_height1,prop_sigover_height1,nsig_height1_5e_4,prop_over_height1_5e_4,prop_sigover_height1_5e_4)

## Height 2
summary_data_height_2 <-  read.table('real_data/data/summary_data_height_2.txt',header=TRUE)

sig_height2 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-8,]
nsig_height2 <- nrow(sig_height2)
prop_over_height2 <- sum(abs(sig_height2$beta)>abs(sig_height2$beta_rep))/nrow(sig_height2)
prop_sigover_height2 <- sum(abs(sig_height2$beta)>(abs(sig_height2$beta_rep)+1.96*sig_height2$se))/nrow(sig_height2)

sig_height2_5e_4 <- summary_data_height_2[2*(pnorm(abs(summary_data_height_2$beta/summary_data_height_2$se), lower.tail=FALSE)) < 5e-4,]
nsig_height2_5e_4 <- nrow(sig_height2_5e_4)
prop_over_height2_5e_4 <- sum(abs(sig_height2_5e_4$beta)>abs(sig_height2_5e_4$beta_rep))/nrow(sig_height2_5e_4)
prop_sigover_height2_5e_4 <- sum(abs(sig_height2_5e_4$beta)>(abs(sig_height2_5e_4$beta_rep)+1.96*sig_height2_5e_4$se))/nrow(sig_height2_5e_4)

values_height2 <- c(nsig_height2,prop_over_height2,prop_sigover_height2,nsig_height2_5e_4,prop_over_height2_5e_4,prop_sigover_height2_5e_4)

################################################################################
## Combine above results and save them in results file

# Table 1: initial exploration
values <- rbind(values_BMI1,values_BMI2, values_T2D1, values_T2D2, values_height1, values_height2)
values <- round(values[,1:6],4)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
values <- cbind(GWAS,values)
colnames(values) <- c("GWAS", "nsig 5e-8", "overest 5e-8","sig_overest 5e-8","nsig 5e-4", "overest 5e-4","sig_overest 5e-4")
rownames(values) <- NULL
write.csv(values,"real_data/results/realdata_Table1.txt", row.names = FALSE)

################################################################################
################################################################################

## PART B)

## 1) BMI 1:
## COMPUTING ESTIMATED MSE AND BIAS FOR BMI 1 AT THRESHOLDS 5e-8 AND 5e-4

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

## BIAS (5e-8 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_1$beta_rep

bias_BMI1_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-8 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_bmi_1$beta_rep

bias_BMI1_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-4 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_1$beta_rep

bias_BMI1_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-4 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]
summary_data_bmi_1 <- dplyr::arrange(summary_data_bmi_1,dplyr::desc(abs(summary_data_bmi_1$beta/summary_data_bmi_1$se)))
summary_data_sig_bmi_1 <- summary_data_bmi_1[summary_data_bmi_1$beta/summary_data_bmi_1$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_bmi_1$beta_rep

bias_BMI1_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

################################################################################

## 2) BMI 2:
## COMPUTING ESTIMATED MSE AND BIAS FOR BMI 2 AT THRESHOLDS 5e-8 AND 5e-4

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

## BIAS (5e-8 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_2$beta_rep

bias_BMI2_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-8 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_bmi_2$beta_rep

bias_BMI2_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-4 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_bmi_2$beta_rep

bias_BMI2_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-4 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]
summary_data_bmi_2 <- dplyr::arrange(summary_data_bmi_2,dplyr::desc(abs(summary_data_bmi_2$beta/summary_data_bmi_2$se)))
summary_data_sig_bmi_2 <- summary_data_bmi_2[summary_data_bmi_2$beta/summary_data_bmi_2$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_bmi_2$beta_rep

bias_BMI2_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

################################################################################

## 3) T2D 1:
## COMPUTING ESTIMATED MSE AND BIAS FOR T2D 1 AT THRESHOLDS 5e-8 AND 5e-4

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

## BIAS (5e-8 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]
summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_1$beta_rep

bias_T2D1_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-8 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]
summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_T2D_1$beta_rep

bias_T2D1_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-4 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]
summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_1$beta_rep

bias_T2D1_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-4 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]
summary_data_T2D_1 <- dplyr::arrange(summary_data_T2D_1,dplyr::desc(abs(summary_data_T2D_1$beta/summary_data_T2D_1$se)))
summary_data_sig_T2D_1 <- summary_data_T2D_1[summary_data_T2D_1$beta/summary_data_T2D_1$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_T2D_1$beta_rep

bias_T2D1_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-8 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]
summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_2$beta_rep

bias_T2D2_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-8 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]
summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_T2D_2$beta_rep

bias_T2D2_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-4 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]
summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_T2D_2$beta_rep

bias_T2D2_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-4 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]
summary_data_T2D_2 <- dplyr::arrange(summary_data_T2D_2,dplyr::desc(abs(summary_data_T2D_2$beta/summary_data_T2D_2$se)))
summary_data_sig_T2D_2 <- summary_data_T2D_2[summary_data_T2D_2$beta/summary_data_T2D_2$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_T2D_2$beta_rep

bias_T2D2_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-8 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_1$beta_rep

bias_height1_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-8 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_height_1$beta_rep

bias_height1_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-4 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_1$beta_rep

bias_height1_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-4 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]
summary_data_height_1 <- dplyr::arrange(summary_data_height_1,dplyr::desc(abs(summary_data_height_1$beta/summary_data_height_1$se)))
summary_data_sig_height_1 <- summary_data_height_1[summary_data_height_1$beta/summary_data_height_1$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_height_1$beta_rep

bias_height1_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-8 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-8/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-8/2, lower.tail=FALSE),]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se > qnorm(5e-8/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_2$beta_rep

bias_height2_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-8 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-8/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-8/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-8/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-8/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-8/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-8/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-8/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-8/2),]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se < qnorm(5e-8/2),]
true_beta <- summary_data_sig_height_2$beta_rep

bias_height2_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

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

## BIAS (5e-4 - positive)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se > qnorm(5e-4/2, lower.tail=FALSE),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se > qnorm(5e-4/2, lower.tail=FALSE),]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se > qnorm(5e-4/2, lower.tail=FALSE),]
true_beta <- summary_data_sig_height_2$beta_rep

bias_height2_5e_4_up <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

## BIAS (5e-4 - negative)
out_CL_sig <- out_CL[out_CL$beta/out_CL$se < qnorm(5e-4/2,),]
out_EB_sig <- out_EB[out_EB$beta/out_EB$se < qnorm(5e-4/2),]
out_EB_df_sig <- out_EB_df[out_EB_df$beta/out_EB_df$se < qnorm(5e-4/2),]
out_EB_scam_sig <- out_EB_scam[out_EB_scam$beta/out_EB_scam$se < qnorm(5e-4/2),]
out_EB_gampo_sig <- out_EB_gam_po[out_EB_gam_po$beta/out_EB_gam_po$se < qnorm(5e-4/2),]
out_EB_gamnb_sig <- out_EB_gam_nb[out_EB_gam_nb$beta/out_EB_gam_nb$se < qnorm(5e-4/2),]
out_FIQT_sig <- out_FIQT[out_FIQT$beta/out_FIQT$se < qnorm(5e-4/2),]
out_BR_ss_sig <- out_BR[out_BR$beta/out_BR$se < qnorm(5e-4/2),]
summary_data_height_2 <- dplyr::arrange(summary_data_height_2,dplyr::desc(abs(summary_data_height_2$beta/summary_data_height_2$se)))
summary_data_sig_height_2 <- summary_data_height_2[summary_data_height_2$beta/summary_data_height_2$se < qnorm(5e-4/2),]
true_beta <- summary_data_sig_height_2$beta_rep

bias_height2_5e_4_down <- data.frame(naive = mean((out_EB_sig$beta - true_beta)), CL1 = mean((out_CL_sig$beta.cl1 - true_beta)), CL2 = mean((out_CL_sig$beta.cl2 - true_beta)), CL3 = mean((out_CL_sig$beta.cl3 - true_beta)), EB = mean((out_EB_sig$beta_EB - true_beta)),  EB_df = mean((out_EB_df_sig$beta_EB - true_beta)), EB_scam = mean((out_EB_scam_sig$beta_EB - true_beta)), EB_gam_po = mean((out_EB_gampo_sig$beta_EB - true_beta)), EB_gam_nb = mean((out_EB_gamnb_sig$beta_EB - true_beta)), boot = mean((out_BR_ss_sig$beta_BR_ss - true_beta)), FIQT = mean((out_FIQT_sig$beta_FIQT - true_beta)))

###############################################################################

## Combine above results and save them in results file

# Table 2: MSE 5e-8
mse_5e_8 <- rbind(mse_BMI1,mse_BMI2, mse_T2D1, mse_T2D2, mse_height1, mse_height2)
mse_5e_8 <- round(mse_5e_8[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_8 <- cbind(GWAS,mse_5e_8)
write.csv(mse_5e_8,"real_data/results/mse_5e_8.txt", row.names = FALSE)

# S12 Table: MSE 5e-4
mse_5e_4 <- rbind(mse_BMI1_5e_4,mse_BMI2_5e_4, mse_T2D1_5e_4, mse_T2D2_5e_4, mse_height1_5e_4, mse_height2_5e_4)
mse_5e_4 <- round(mse_5e_4[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
mse_5e_4 <- cbind(GWAS,mse_5e_4)
write.csv(mse_5e_4,"real_data/results/mse_5e_4.txt", row.names = FALSE)

# S10 Table: BIAS 5e-8 (positive)
bias_5e_8_up <- rbind(bias_BMI1_up,bias_BMI2_up, bias_T2D1_up, bias_T2D2_up, bias_height1_up, bias_height2_up)
bias_5e_8_up <- round(bias_5e_8_up[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_8_up <- cbind(GWAS,bias_5e_8_up)
write.csv(bias_5e_8_up,"real_data/results/bias_5e_8_positive.txt", row.names = FALSE)

# S11 Table: BIAS 5e-8 (negative)
bias_5e_8_down <- rbind(bias_BMI1_down,bias_BMI2_down, bias_T2D1_down, bias_T2D2_down, bias_height1_down, bias_height2_down)
bias_5e_8_down <- round(bias_5e_8_down[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_8_down <- cbind(GWAS,bias_5e_8_down)
write.csv(bias_5e_8_down,"real_data/results/bias_5e_8_negative.txt", row.names = FALSE)

# S13 Table: BIAS 5e-4 (positive)
bias_5e_4_up <- rbind(bias_BMI1_5e_4_up,bias_BMI2_5e_4_up, bias_T2D1_5e_4_up, bias_T2D2_5e_4_up, bias_height1_5e_4_up, bias_height2_5e_4_up)
bias_5e_4_up <- round(bias_5e_4_up[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_4_up <- cbind(GWAS,bias_5e_4_up)
write.csv(bias_5e_4_up,"real_data/results/bias_5e_4_positive.txt", row.names = FALSE)

# S14 Table: BIAS 5e-4 (negative)
bias_5e_4_down <- rbind(bias_BMI1_5e_4_down,bias_BMI2_5e_4_down, bias_T2D1_5e_4_down, bias_T2D2_5e_4_down, bias_height1_5e_4_down, bias_height2_5e_4_down)
bias_5e_4_down <- round(bias_5e_4_down[,1:11],5)
GWAS <- c("BMI 1", "BMI 2", "T2D 1", "T2D 2", "Height 1", "Height 2")
bias_5e_4_down <- cbind(GWAS,bias_5e_4_down)
write.csv(bias_5e_4_down,"real_data/results/bias_5e_4_negative.txt", row.names = FALSE)

###############################################################################

