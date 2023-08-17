## WINNER'S CURSE SIMULATION STUDY SCRIPT 6 - EVALUATING METHODS (IND2)::

## This script evaluates and compares several winner's curse correction methods
## using simulated sets of summary statistics, in which SNPs are INDEPENDENT.
## The simulated data sets correspond to three different settings: a
## quantitative trait with bimodal effect size distribution, a quantitative
## trait with skewed effect size distribution, and a binary trait with normal
## effect size distribution. At two significance thresholds, 5e-8 and 5e-4, this
## script compares the following methods using the following evaluation metrics:

## Methods evaluated: 
## 1. original empirical Bayes method with tail restriction
## 2. FDR Inverse Quantile Transformation
## 3. proposed bootstrap approach
## 4. conditional likelihood estimators of Ghosh et al. (2008)  

## Evaluation metrics: 
## 1. flb - fraction of significant SNPs less biased 
## 2. mse - improvement in average MSE for significant SNPs 
## 3. rmse - improvement in average RMSE for significant SNPs
## 4. bias_up - average bias of significant SNPs with positive estimated effects
## 5. bias_down - average bias of significant SNPs with negative estimated effects
## 6. rel_mse - relative improvement in average MSE for significant SNPs

## Outputs: 
## 1. bim_5e-8_100sim_ave.csv
## 2. bim_5e-8_100sim_all.csv
## 3. bim_5e-4_100sim_ave.csv
## 4. bim_5e-4_100sim_all.csv
## 5. skew_5e-8_100sim_ave.csv
## 6. skew_5e-8_100sim_all.csv
## 7. skew_5e-4_100sim_ave.csv
## 8. skew_5e-4_100sim_all.csv
## 9. bin_5e-8_100sim_ave.csv
## 10. bin_5e-8_100sim_all.csv
## 11. bin_5e-4_100sim_ave.csv
## 12. bin_5e-4_100sim_all.csv

################################################################################

## Total number of simulations:
tot_sim <- 100
## Fixed total number of SNPs:
n_snps <- 10^6

## Set of scenarios to be tested
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  n_samples = c(30000,300000),
  h2 = c(0.3,0.8),
  prop_effect = c(0.01,0.001),
  S = c(-1, 0, 1)
)

################################################################################
## PART 1A) QUANTITATIVE TRAIT, BIMODAL DISTRIBUTION, THRESHOLD 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bim(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  
  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-8)
  
  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  
  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-8)
  
  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-8)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-8)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-8)
  
  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=5)
  
  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)
  
  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-8)
  
  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=2)
  mse_naive <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=2)
  bias_naive_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=2)
  bias_naive_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=2)
  
  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################
## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## boot:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)


## cl1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## cl2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## cl3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## rep:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)), rep("naive", (nrow(sim_params)/tot_sim)))
write.csv(results_all,"simulations/results/bim_5e-8_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep, res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)), rep("naive",nrow(sim_params)))
write.csv(results_all,"simulations/results/bim_5e-8_100sim_all.csv")

print("PART 1A) complete!")

################################################################################
## PART 1B) QUANTITATIVE TRAIT, BIMODAL DISTRIBUTION, THRESHOLD 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bim(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  
  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-4)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-4)
  
  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  
  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-4)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-4)
  
  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-4)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-4)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-4)
  
  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=5)
  
  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=6)
  
  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-4)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-4)
  
  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=2)
  mse_naive <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=2)
  bias_naive_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=2)
  bias_naive_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=2)
  
  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################
## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## boot:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)


## cl1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## cl2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## cl3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## rep:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)), rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"simulations/results/bim_5e-4_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep, res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)), rep("naive",nrow(sim_params)))
write.csv(results_all,"simulations/results/bim_5e-4_100sim_all.csv")

print("PART 1B) complete!")

################################################################################
## PART 2A) QUANTITATIVE TRAIT, SKEWED DISTRIBUTION, THRESHOLD 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_exp(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  
  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-8)
  
  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  
  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-8)
  
  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-8)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-8)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-8)
  
  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=5)
  
  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)
  
  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-8)
  
  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=2)
  mse_naive <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=2)
  bias_naive_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=2)
  bias_naive_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=2)
  
  
  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################
## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## boot:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)


## cl1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## cl2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## cl3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## rep:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)), rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/skew_5e-8_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep, res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)), rep("naive",nrow(sim_params)))
write.csv(results_all,"results/skew_5e-8_100sim_all.csv")

print("PART 2A) complete!")

################################################################################
## PART 2B) QUANTITATIVE TRAIT, SKEWED DISTRIBUTION, THRESHOLD 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_exp(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  
  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-4)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-4)
  
  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  
  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-4)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-4)
  
  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-4)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-4)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-4)
  
  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=5)
  
  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=6)
  
  ## replication
  ss2 <- data.frame(true_beta=ss$true_beta,se=ss$rep_se)
  rep_stats <- simulate_est(ss2)
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-4)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-4)
  
  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=2)
  mse_naive <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=2)
  bias_naive_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=2)
  bias_naive_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=2)
  
  
  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################
## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## boot:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## cl1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## cl2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## cl3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## rep:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)), rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"simulations/results/skew_5e-4_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep, res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)), rep("naive",nrow(sim_params)))
write.csv(results_all,"simulations/results/skew_5e-4_100sim_all.csv")

################################################################################
## PART 3A) BINARY TRAIT, NORMAL DISTRIBUTION, THRESHOLD 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bin(H2=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss)
  disc_stats <- disc_stats[!is.na(disc_stats$beta) & !is.na(rep_stats$beta),]
  rep_stats <- rep_stats[!is.na(disc_stats$beta) & !is.na(rep_stats$beta),]
  
  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-8)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-8)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-8)
  
  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-8)
  
  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-8)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-8)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-8)
  
  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-8)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-8)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-8)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-8)
  
  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=5)
  
  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=6)
  
  ## replication
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-8)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-8)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-8)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-8)
  
  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-8,i=2)
  mse_naive <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-8,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-8,i=2)
  bias_naive_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-8,i=2)
  bias_naive_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-8,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-8,i=2)
  
  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################
## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## boot:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)


## cl1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## cl2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## cl3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## rep:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)), rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"results/bin_5e-8_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep, res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)), rep("naive",nrow(sim_params)))
write.csv(results_all,"results/bin_5e-8_100sim_all.csv")

print("PART 3A) complete!")

################################################################################
## PART 3B) BINARY TRAIT, NORMAL DISTRIBUTION, THRESHOLD 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  ss <- simulate_ss_bin(H2=h2,Pi=prop_effect,nid=n_samples,sc=S)
  disc_stats <- simulate_est(ss)
  rep_stats <- simulate_est(ss)
  disc_stats <- disc_stats[!is.na(disc_stats$beta) & !is.na(rep_stats$beta),]
  rep_stats <- rep_stats[!is.na(disc_stats$beta) & !is.na(rep_stats$beta),]
  
  ## Empirical Bayes:
  out_EB <- empirical_bayes(disc_stats)
  flb_EB <- frac_sig_less_bias(out_EB,ss$true_beta,i=4,alpha=5e-4)
  mse_EB <- mse_sig_improve(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rmse_EB <- mse_sig_improve_root(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_up <- bias_sig_up(out_EB,ss$true_beta,i=4,alpha=5e-4)
  bias_EB_down <- bias_sig_down(out_EB,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_EB <- mse_sig_improve_per(out_EB,ss$true_beta,i=4,alpha=5e-4)
  
  ## FIQT:
  out_FIQT <- FDR_IQT(disc_stats)
  flb_FIQT <- frac_sig_less_bias(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  mse_FIQT <- mse_sig_improve(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rmse_FIQT <- mse_sig_improve_root(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_up <- bias_sig_up(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  bias_FIQT_down <- bias_sig_down(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_FIQT <- mse_sig_improve_per(out_FIQT,ss$true_beta,i=4,alpha=5e-4)
  
  ## Bootstrap:
  out_BR <- BR_ss(disc_stats)
  flb_BR <- frac_sig_less_bias(out_BR,ss$true_beta,i=4,alpha=5e-4)
  mse_BR <- mse_sig_improve(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rmse_BR <- mse_sig_improve_root(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_up <- bias_sig_up(out_BR,ss$true_beta,i=4,alpha=5e-4)
  bias_BR_down <- bias_sig_down(out_BR,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_BR <- mse_sig_improve_per(out_BR,ss$true_beta,i=4,alpha=5e-4)
  
  ## cl1:
  out_cl <- conditional_likelihood_ind(disc_stats,alpha=5e-4)
  flb_cl1 <- frac_sig_less_bias(out_cl,ss$true_beta,i=4,alpha=5e-4)
  mse_cl1 <- mse_sig_improve(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rmse_cl1 <- mse_sig_improve_root(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_up <- bias_sig_up(out_cl,ss$true_beta,i=4,alpha=5e-4)
  bias_cl1_down <- bias_sig_down(out_cl,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_cl1 <- mse_sig_improve_per(out_cl,ss$true_beta,i=4,alpha=5e-4)
  
  ## cl2:
  flb_cl2 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=5)
  mse_cl2 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rmse_cl2 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=5)
  bias_cl2_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=5)
  rel_mse_cl2 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=5)
  
  ## cl3:
  flb_cl3 <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=6)
  mse_cl3 <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rmse_cl3 <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=6)
  bias_cl3_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=6)
  rel_mse_cl3 <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=6)
  
  ## replication
  out_rep <- cbind(disc_stats,rep_stats$beta)
  flb_rep <- frac_sig_less_bias(out_rep,ss$true_beta,i=4,alpha=5e-4)
  mse_rep <- mse_sig_improve(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rmse_rep <- mse_sig_improve_root(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_up <- bias_sig_up(out_rep,ss$true_beta,i=4,alpha=5e-4)
  bias_rep_down <- bias_sig_down(out_rep,ss$true_beta,i=4,alpha=5e-4)
  rel_mse_rep <- mse_sig_improve_per(out_rep,ss$true_beta,i=4,alpha=5e-4)
  
  ## naive - only interested in bias (improvement should be 0):
  flb_naive <- frac_sig_less_bias(out_cl,ss$true_beta,alpha=5e-4,i=2)
  mse_naive <- mse_sig_improve(out_cl,ss$true_beta,alpha=5e-4,i=2)
  rmse_naive <- mse_sig_improve_root(out_cl,ss$true_beta,alpha=5e-4,i=2)
  bias_naive_up <- bias_sig_up(out_cl,ss$true_beta,alpha=5e-4,i=2)
  bias_naive_down <- bias_sig_down(out_cl,ss$true_beta,alpha=5e-4,i=2)
  rel_mse_naive <- mse_sig_improve_per(out_cl,ss$true_beta,alpha=5e-4,i=2)
  
  
  return(c(flb_EB,mse_EB,rmse_EB,bias_EB_up,bias_EB_down,rel_mse_EB,flb_FIQT,mse_FIQT,rmse_FIQT,bias_FIQT_up,bias_FIQT_down,rel_mse_FIQT,flb_BR,mse_BR,rmse_BR,bias_BR_up,bias_BR_down,rel_mse_BR,flb_cl1,mse_cl1,rmse_cl1,bias_cl1_up,bias_cl1_down,rel_mse_cl1,flb_cl2,mse_cl2,rmse_cl2,bias_cl2_up,bias_cl2_down,rel_mse_cl2,flb_cl3,mse_cl3,rmse_cl3,bias_cl3_up,bias_cl3_down,rel_mse_cl3,flb_rep,mse_rep,rmse_rep,bias_rep_up,bias_rep_down,rel_mse_rep,flb_naive,mse_naive,rmse_naive,bias_naive_up,bias_naive_down,rel_mse_naive))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

################################################################################
## Organising results:

## Empirical Bayes:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][1]
  mse[i] <- res[[i]][2]
  rmse[i] <- res[[i]][3]
  bias_up[i] <- res[[i]][4]
  bias_down[i] <- res[[i]][5]
  rel_mse[i] <- res[[i]][6]
}
res_EB <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_EB <- ave_results(res_EB,tot_sim)

## FIQT:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][7]
  mse[i] <- res[[i]][8]
  rmse[i] <- res[[i]][9]
  bias_up[i] <- res[[i]][10]
  bias_down[i] <- res[[i]][11]
  rel_mse[i] <- res[[i]][12]
}
res_FIQT <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_FIQT <- ave_results(res_FIQT,tot_sim)

## boot:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][13]
  mse[i] <- res[[i]][14]
  rmse[i] <- res[[i]][15]
  bias_up[i] <- res[[i]][16]
  bias_down[i] <- res[[i]][17]
  rel_mse[i] <- res[[i]][18]
}
res_BR <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_BR <- ave_results(res_BR,tot_sim)

## cl1:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][19]
  mse[i] <- res[[i]][20]
  rmse[i] <- res[[i]][21]
  bias_up[i] <- res[[i]][22]
  bias_down[i] <- res[[i]][23]
  rel_mse[i] <- res[[i]][24]
}
res_cl1 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl1 <- ave_results(res_cl1,tot_sim)

## cl2:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][25]
  mse[i] <- res[[i]][26]
  rmse[i] <- res[[i]][27]
  bias_up[i] <- res[[i]][28]
  bias_down[i] <- res[[i]][29]
  rel_mse[i] <- res[[i]][30]
}
res_cl2 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl2 <- ave_results(res_cl2,tot_sim)

## cl3:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][31]
  mse[i] <- res[[i]][32]
  rmse[i] <- res[[i]][33]
  bias_up[i] <- res[[i]][34]
  bias_down[i] <- res[[i]][35]
  rel_mse[i] <- res[[i]][36]
}
res_cl3 <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_cl3 <- ave_results(res_cl3,tot_sim)

## rep:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][37]
  mse[i] <- res[[i]][38]
  rmse[i] <- res[[i]][39]
  bias_up[i] <- res[[i]][40]
  bias_down[i] <- res[[i]][41]
  rel_mse[i] <- res[[i]][42]
}
res_rep <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_rep <- ave_results(res_rep,tot_sim)

## naive:
flb <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
rmse <- c(rep(0,nrow(sim_params)))
bias_up <- c(rep(0,nrow(sim_params)))
bias_down <- c(rep(0,nrow(sim_params)))
rel_mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  flb[i] <- res[[i]][43]
  mse[i] <- res[[i]][44]
  rmse[i] <- res[[i]][45]
  bias_up[i] <- res[[i]][46]
  bias_down[i] <- res[[i]][47]
  rel_mse[i] <- res[[i]][48]
}
res_naive <- cbind(sim_params,flb,mse,rmse,bias_up,bias_down,rel_mse)
ave_res_naive <- ave_results(res_naive,tot_sim)

## Combine all results:
results_all <- rbind(ave_res_EB,ave_res_FIQT,ave_res_BR,ave_res_cl1,ave_res_cl2,ave_res_cl3, ave_res_rep, ave_res_naive)
results_all$method <- c(rep("EB",(nrow(sim_params)/tot_sim)),rep("FIQT",(nrow(sim_params)/tot_sim)),rep("BR",(nrow(sim_params)/tot_sim)),rep("cl1",(nrow(sim_params)/tot_sim)),rep("cl2",(nrow(sim_params)/tot_sim)),rep("cl3",(nrow(sim_params)/tot_sim)), rep("rep",(nrow(sim_params)/tot_sim)), rep("naive",(nrow(sim_params)/tot_sim)))
write.csv(results_all,"simulations/results/bin_5e-4_100sim_ave.csv")

## Combine all results:
results_all <- rbind(res_EB,res_FIQT,res_BR,res_cl1,res_cl2,res_cl3,res_rep, res_naive)
results_all$method <- c(rep("EB",nrow(sim_params)),rep("FIQT",nrow(sim_params)),rep("BR",nrow(sim_params)),rep("cl1",nrow(sim_params)),rep("cl2",nrow(sim_params)),rep("cl3",nrow(sim_params)),rep("rep",nrow(sim_params)), rep("naive",nrow(sim_params)))
write.csv(results_all,"simulations/results/bin_5e-4_100sim_all.csv")

print("PART 3B) complete!")

################################################################################

