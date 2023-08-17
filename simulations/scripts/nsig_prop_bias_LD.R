## WINNER'S CURSE SIMULATION STUDY SCRIPT 3 - INITIAL EXPLORATION (LD):

## This script performs a preliminary investigation of the simulated sets of
## summary statistics, in which a simple correlation structure has been imposed
## on independent blocks of 100 SNPs. The simulated data sets also correspond to
## a quantitative trait with a normal effect size distribution. At two
## significance thresholds, 5e-8 and 5e-4, this script obtains the number of
## significant SNPs, the proportion of these SNPs for which their association
## estimate is more extreme than their true effect size, the proportion of these
## SNPs which are significantly overestimated and the mean square error (MSE) of
## significant SNPs for each simulated data set.

## Outputs: 
## 1. nsig_prop_bias_5e_8_LD_all.csv
## 2. nsig_prop_bias_5e-8_LD.csv (used to obtain S1 Table)
## 3. nsig_prop_bias_5e_4_LD_all.csv
## 4. nsig_prop_bias_5e-4_LD.csv

################################################################################

## Total number of simulations: 100
tot_sim <- 100
## Fixed total number of SNPs:
n_snps <- 10^6

## Set of scenarios to be tested
sim_params <- expand.grid(
  sim = c(1:tot_sim),
  n_samples = c(30000,300000),
  h2 = c(0.3,0.8),
  prop_effect = c(0.01, 0.001),
  S = c(0)
)

################################################################################
## PART 1A) THRESHOLD 5e-8

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  out <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  snp_sig <- out[abs(out$beta/out$se) > qnorm((5e-8)/2, lower.tail=FALSE),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
  }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(out$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(out$true_beta[snp_sig$rsid]) + 1.96*out$se[snp_sig$rsid]))/n_sig
    mse <- mean((out$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

n_sig <- c(rep(0,nrow(sim_params)))
prop_bias <- c(rep(0,nrow(sim_params)))
prop_x <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  n_sig[i] <- res[[i]][1]
  prop_bias[i] <- res[[i]][2]
  prop_x[i] <- res[[i]][3]
  mse[i] <- res[[i]][4]
}

results <- cbind(sim_params,n_sig,prop_bias,prop_x,mse)
ave_res <- ave_results1(results,tot_sim)

write.csv(results, "simulations/results/nsig_prop_bias_5e_8_LD_all.csv")
write.csv(ave_res, "simulations/results/nsig_prop_bias_5e-8_LD.csv")

print("PART 1A) complete!")

################################################################################
## PART 1B) THRESHOLD 5e-4

set.seed(1998)

run_sim <- function(n_samples, h2, prop_effect, S,sim)
{
  out <- simulate_ss_ld(H=h2,Pi=prop_effect,nid=n_samples,sc=S)
  snp_sig <- out[abs(out$beta/out$se) > qnorm((5e-4)/2, lower.tail=FALSE),]
  n_sig <- nrow(snp_sig)
  if (n_sig == 0){
    prop_bias <- -1
    mse <- -1
    prop_x <- -1
  }else{
    prop_bias <- sum(abs(snp_sig$beta) > abs(out$true_beta[snp_sig$rsid]))/n_sig
    prop_x <- sum(abs(snp_sig$beta) > (abs(out$true_beta[snp_sig$rsid]) + 1.96*out$se[snp_sig$rsid]))/n_sig
    mse <- mean((out$true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  }
  return(c(n_sig,prop_bias,prop_x,mse))
}
res <- mclapply(1:nrow(sim_params), function(i){
  #print(paste(round(i*100/nrow(sim_params), 2),"%"))
  do.call(run_sim, args=as.list(sim_params[i,]))}, mc.cores=1)

n_sig <- c(rep(0,nrow(sim_params)))
prop_bias <- c(rep(0,nrow(sim_params)))
prop_x <- c(rep(0,nrow(sim_params)))
mse <- c(rep(0,nrow(sim_params)))
for (i in 1:nrow(sim_params)){
  n_sig[i] <- res[[i]][1]
  prop_bias[i] <- res[[i]][2]
  prop_x[i] <- res[[i]][3]
  mse[i] <- res[[i]][4]
}

results <- cbind(sim_params,n_sig,prop_bias,prop_x,mse)
ave_res <- ave_results1(results,tot_sim)

write.csv(results, "simulations/results/nsig_prop_bias_5e_4_LD_all.csv")
write.csv(ave_res, "simulations/results/nsig_prop_bias_5e-4_LD.csv")

print("PART 1B) complete!")

################################################################################
