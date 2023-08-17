## WINNER'S CURSE SIMULATION STUDY SCRIPT 2 - FUNCTIONS:

## This script provides all functions required for simulations. 

################################################################################
# 1. simulate_ss_ld() - Given values for heritability (H), polygenicity (Pi),
# sample size (nid), selection coefficient (S) and given R, the 100x100 LD matrix
# of inter-genotype correlations, 'true' values of effect size, estimated effect
# sizes and corresponding standard errors are simulated for a quantitative trait 
# with a simple correlation structure.

# As R is fixed, it is defined outside the function as follows.

s <- function(n,x=0.9825){
  vector <- c()
  for (i in 1:n){vector[i] <- x^(i-1)}
  return(vector)
}

vec <- function(tot){
  vector <- c()
  for (i in 1:tot){vector <- c(vector,s(100-(i-1)))}
  return(vector)
}

m1 <- matrix(NA, 100, 100)
m1[lower.tri(m1, diag=TRUE)] <- vec(100)
m2 <- t(m1)
m2[lower.tri(m2, diag=TRUE)] <- vec(100)
R <- m2
R_sqrt <- sqrtm(R) # remains constant

simulate_ss_ld <- function(H,Pi,nid,sc,cormat=R,cormatsq=R_sqrt){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf[1:effect_snps]*(1-maf[1:effect_snps]))^sc))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta_all <- c(rep(0,n_snps))
  maf_re <- c(rep(0,n_snps))
  
  positions <- sample(1:n_snps, effect_snps, replace=F)
  for (i in 1:effect_snps){
    true_beta_all[positions[i]] <- true_beta[i]
    maf_re[positions[i]] <- maf[i]
  }
  
  vec <- seq(1,10^6)
  vec_new <- vec[! vec %in% positions]
  for (i in (effect_snps+1):n_snps){
    maf_re[vec_new[i-effect_snps]] <- maf[i]
  }
  
  mu_all <- numeric(n_snps)
  beta_hat_all <- numeric(n_snps)
  se_all <- numeric(n_snps)
  
  for (i in 1:10000){
    D_1 <- diag(x = 1/(sqrt(nid*2*maf_re[((i-1)*100+1):(i*100)]*(1-maf_re[((i-1)*100+1):(i*100)]))), nrow=100, ncol=100)
    D_2 <- diag(x = (sqrt(nid*2*maf_re[((i-1)*100+1):(i*100)]*(1-maf_re[((i-1)*100+1):(i*100)]))), nrow=100, ncol=100)
    A <- D_1 %*% cormat %*% D_2
    B <- cormatsq %*% D_1
    mu_all[((i-1)*100+1):(i*100)] <- A %*% true_beta_all[((i-1)*100+1):(i*100)]
    beta_hat_all[((i-1)*100+1):(i*100)] <- mu_all[((i-1)*100+1):(i*100)] + ( B %*% rnorm(100))
    se_all[((i-1)*100+1):(i*100)] <- diag(D_1)
  }
  
  summary_stats <- data.frame(rsid=seq(1,n_snps),beta=beta_hat_all,se=se_all,true_beta=mu_all)
  return(summary_stats)
}

################################################################################
# 2. ave_results1() - Takes the 'long' results data frame which contains number
# of significant SNPs and proportion of these SNP that are 'biased' for each
# simulation and outputs means and standard deviations for both quantities
# over all simulations of the same parameter values.

# NOTE: Statistics for proportion of biased significant SNPs are only obtained
# over simulations in which at least one significant SNP was detected.

ave_results1 <- function(res_vec, n_sim){
  ave_nsig <- c(rep(0,nrow(res_vec)/n_sim))
  sd_nsig <- c(rep(0,nrow(res_vec)/n_sim))
  ave_pb <- c(rep(0,nrow(res_vec)/n_sim))
  sd_pb <- c(rep(0,nrow(res_vec)/n_sim))
  ave_px <- c(rep(0,nrow(res_vec)/n_sim))
  sd_px <- c(rep(0,nrow(res_vec)/n_sim))
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_nsig[i] <- mean(res_vec$n_sig[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_nsig[i] <- sd(res_vec$n_sig[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_pb[i] <- mean(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    sd_pb[i] <- sd(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_bias[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    ave_px[i] <- mean(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    sd_px[i] <- sd(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$prop_x[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    ave_mse[i] <- mean(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
    sd_mse[i] <- sd(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)][which(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)] != -1)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$n_sig <- ave_nsig
  res_vec_ave$n_sig_sd <- sd_nsig
  res_vec_ave$prop_bias <- ave_pb
  res_vec_ave$prop_bias_sd <- sd_pb
  res_vec_ave$prop_x <- ave_px
  res_vec_ave$prop_x_sd <- sd_px
  res_vec_ave$mse <- ave_mse
  res_vec_ave$mse_sd <- sd_mse
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}

################################################################################
# 3. frac_sig_less_bias() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the fraction
# of the significant SNPs that are now less biased, i.e. closer to their true 
# value, due to method application is computed. 

frac_sig_less_bias <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  flb <- sum(abs(true_beta[snp_sig$rsid] - snp_sig$beta)>abs(true_beta[snp_sig$rsid] - snp_sig[,i]))/length(snp_sig$rsid)
  return(flb)
}

################################################################################
# 4. mse_sig_improve() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the 
# improvement in average MSE for significant SNPs due to method application is
# computed.

mse_sig_improve <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)-mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)
  return(mse_sig_improve)
}

################################################################################
# 5. mse_sig_improve_root() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the 
# improvement in average RMSE for significant SNPs due to method application is
# computed.

mse_sig_improve_root <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- sqrt(mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2))-sqrt(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))
  return(mse_sig_improve)
}

################################################################################
# 6. bias_sig_up() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the 
# average bias of significant SNPs, with positive estimated effect sizes, is
# computed.

bias_sig_up <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  bias_sig <- mean((snp_sig[,i]-true_beta[snp_sig$rsid]))
  return(bias_sig)
}

################################################################################
# 7. bias_sig_down() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the 
# average bias of significant SNPs, with negative estimated effect sizes, is
# computed.

bias_sig_down <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[(out$beta/out$se) < qnorm((alpha)/2),]
  if (nrow(snp_sig) == 0){return(100)}
  bias_sig <- mean((snp_sig[,i]-true_beta[snp_sig$rsid]))
  return(bias_sig)
}

################################################################################
# 8. rmse_sig_improve_per() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the relative
# improvement in average RMSE for significant SNPs due to method application is
# computed.

rmse_sig_improve_per <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- (sqrt(mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2))-sqrt(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)))/(sqrt(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2)))
  return(mse_sig_improve)
}

################################################################################
# 9. conditional_likelihood_old() - Older version of the function 
# conditional_likelihood() from the 'winnerscurse' package which is better 
# suited for simulations as it doesn't include warning message. Suitable for use
# with LD simulations. 

conditional_likelihood_old <- function(summary_data, alpha=5e-8){
  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))
  
  z <- summary_data$beta/summary_data$se
  p_val <- 2*(stats::pnorm(abs(z), lower.tail=FALSE))
  summary_data <- cbind(summary_data, z, p_val)
  
  if(sum(summary_data$p_val<alpha) == 0){
    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$z)))
    return(summary_data[,1:4])
  }
  summary_data_sig <- summary_data[summary_data$p_val<alpha,]
  c <- stats::qnorm((alpha)/2, lower.tail=FALSE)
  beta.cl1 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl2 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl3 <- c(rep(0,nrow(summary_data_sig)))
  
  for (i in 1:nrow(summary_data_sig)){
    cond.like <- function (mu) { return ((stats::dnorm(summary_data_sig$z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}
    beta.cl1[i] <- (stats::optimize(cond.like, c(0,summary_data_sig$z[i]), maximum=TRUE)$maximum)*summary_data_sig$se[i]
    
    if(abs(summary_data_sig$z[i]) < 37){
      beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37)$value)/(stats::integrate(cond.like,-37,37)$value))*(summary_data_sig$se[i])
    }else{
      if(abs(summary_data_sig$z[i]) > 100){
        beta.cl2[i] <- summary_data_sig$beta[i]
      }else{
        beta.cl2[i] <- ((stats::integrate(mean.cond.like,-100,100)$value)/(stats::integrate(cond.like,-100,100)$value))*(summary_data_sig$se[i])
      }
    }
    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }
  
  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)
  if(sum(abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) >= 1){
    bad.cl2 <- which((abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) == TRUE)
    bad.cl2.df <- data.frame(rsid = summary_data_sig$rsid[bad.cl2], z = summary_data_sig$z[bad.cl2], z.cl2 = beta.cl2[bad.cl2]/summary_data_sig$se[bad.cl2], beta = summary_data_sig$beta[bad.cl2], beta.cl2 = beta.cl2[bad.cl2], se = summary_data_sig$se[bad.cl2])
    good.cl2.df <- data.frame(rsid = summary_data_sig$rsid[-(bad.cl2)], z = summary_data_sig$z[-(bad.cl2)], z.cl2 = beta.cl2[-(bad.cl2)]/summary_data_sig$se[-(bad.cl2)], beta = summary_data_sig$beta[-bad.cl2], beta.cl2 = beta.cl2[-(bad.cl2)], se = summary_data_sig$se[-(bad.cl2)])
    for(i in 1:nrow(bad.cl2.df)){
      if(length(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]]) == 0 || length(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]]) == 0){bad.cl2.df$z.cl2[i] <- good.cl2.df$z.cl2[which.min(abs(good.cl2.df$z-bad.cl2.df$z[i]))]
      bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]}else{
        min_greater <- min(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]])
        min_smaller <- min(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]])
        bad.cl2.df$z.cl2[i] <- ((abs(min_greater - bad.cl2.df$z[i]))*(good.cl2.df$z.cl2[good.cl2.df$z == min_smaller]) + (abs(bad.cl2.df$z[i] -     min_smaller))*(good.cl2.df$z.cl2[good.cl2.df$z == min_greater]))/(abs(min_greater - min_smaller))
        bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]
      }
      
      summary_data_sig$beta.cl2[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- bad.cl2.df$beta.cl2[i]
      summary_data_sig$beta.cl3[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- (summary_data_sig$beta.cl1[summary_data_sig$rsid == bad.cl2.df$rsid[i]] +bad.cl2.df$beta.cl2[i])/2
    }
  }
  
  summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$z)))
  summary_data_sig <- cbind(summary_data_sig[,1:4],summary_data_sig[,7:9])
  return(summary_data_sig)
}

################################################################################
# 10. average() - Required for use in ave_results(). Obtains mean of evaluation
# metric over only those simulations in which at least one significant SNP has
# been found.

average <- function(ave_vector){
  ave_vector <- ave_vector[which(ave_vector != 100)]
  ave <- mean(ave_vector)
  return(ave)
}

################################################################################
# 11. std() - Required for use in ave_results(). Obtains standard deviation of
# evaluation metric over only those simulations in which at least one significant
# SNP has been found.

std <- function(ave_vector){
  ave_vector <- ave_vector[which(ave_vector != 100)]
  error <- sd(ave_vector)
  return(error)
}

################################################################################
# 12. ave_results() - Takes the 'long' results data frame which contains
# evaluation metrics for each simulation and outputs means and standard
# deviations of these metrics over simulations of the same parameter values, in
# which at least one significant SNP has been found.

ave_results <- function(res_vec, n_sim){
  ave_flb <- c(rep(0,nrow(res_vec)/n_sim))
  sd_flb <- c(rep(0,nrow(res_vec)/n_sim))
  ave_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_mse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_rmse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_rmse <- c(rep(0,nrow(res_vec)/n_sim))
  ave_bias_up <- c(rep(0,nrow(res_vec)/n_sim))
  sd_bias_up <- c(rep(0,nrow(res_vec)/n_sim))
  ave_bias_down <- c(rep(0,nrow(res_vec)/n_sim))
  sd_bias_down <- c(rep(0,nrow(res_vec)/n_sim))
  ave_rel_mse <- c(rep(0,nrow(res_vec)/n_sim))
  sd_rel_mse <- c(rep(0,nrow(res_vec)/n_sim))
  for (i in 1:(nrow(res_vec)/n_sim)){
    ave_flb[i] <- average(res_vec$flb[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_flb[i] <- std(res_vec$flb[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_mse[i] <- average(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_mse[i] <- std(res_vec$mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_rmse[i] <- average(res_vec$rmse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_rmse[i] <- std(res_vec$rmse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_bias_up[i] <- average(res_vec$bias_up[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_bias_up[i] <- std(res_vec$bias_up[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_bias_down[i] <- average(res_vec$bias_down[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_bias_down[i] <- std(res_vec$bias_down[(i*n_sim-(n_sim-1)):(i*n_sim)])
    ave_rel_mse[i] <- average(res_vec$rel_mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
    sd_rel_mse[i] <- std(res_vec$rel_mse[(i*n_sim-(n_sim-1)):(i*n_sim)])
  }
  keep <- seq(from=1,to=nrow(res_vec)-n_sim+1,by=n_sim)
  res_vec_ave <- res_vec[rownames(res_vec) %in% keep,]
  res_vec_ave$flb <- ave_flb
  res_vec_ave$flb_error <- sd_flb
  res_vec_ave$mse <- ave_mse
  res_vec_ave$mse_error <- sd_mse
  res_vec_ave$rmse <- ave_rmse
  res_vec_ave$rmse_error <- sd_rmse
  res_vec_ave$bias_up <- ave_bias_up
  res_vec_ave$bias_up_error <- sd_bias_up
  res_vec_ave$bias_down <- ave_bias_down
  res_vec_ave$bias_down_error <- sd_bias_down
  res_vec_ave$rel_mse <- ave_rel_mse
  res_vec_ave$rel_mse_error <- sd_rel_mse
  res_vec_ave <- subset(res_vec_ave, select = -sim )
  return(res_vec_ave)
}

################################################################################
# 13. simulate_ss() - Given values for heritability (H), polygenicity (Pi),
# sample size (nid) and selection coefficient (S), 'true' values of effect size
# and corresponding standard error are simulated in which effect sizes are
# assumed to be independent and normally distributed.

simulate_ss <- function(H,Pi,nid,sc,rep_nid=1){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf*(1-maf))^sc))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  rep_se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid*nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se,rep_se)
  return(stats)
}


################################################################################
# 14. simulate_est() - Given the output from a function such as simulate_ss(), 
# GWAS summary statistics are simulated. The resulting data frame contains SNP
# ID numbers, estimated effect sizes and standard errors.

simulate_est <- function(stats){
  est <- data.frame(rsid=seq(1,n_snps),beta=rnorm(n=n_snps,mean=stats$true_beta,sd=stats$se),se=stats$se)
  return(est)
}

################################################################################
# 15. conditional_likelihood_ind() - Older version of the function 
# conditional_likelihood() from the 'winnerscurse' package which is better 
# suited for simulations as it doesn't include warning message. Similar to 
# function 9. but suitable for use with independent simulations.

conditional_likelihood_ind <- function(summary_data, alpha=5e-8){
  stopifnot(all(c("rsid", "beta","se") %in% names(summary_data)))
  stopifnot(!all(is.na(summary_data$rsid)) && !all(is.na(summary_data$beta)) && !all(is.na(summary_data$se)))
  stopifnot(is.numeric(summary_data$beta) && is.numeric(summary_data$se))
  stopifnot(!any(duplicated(summary_data$rsid)))
  
  z <- summary_data$beta/summary_data$se
  p_val <- 2*(stats::pnorm(abs(z), lower.tail=FALSE))
  summary_data <- cbind(summary_data, z, p_val)
  
  if(sum(summary_data$p_val<alpha) == 0){
    summary_data <- dplyr::arrange(summary_data,dplyr::desc(abs(summary_data$z)))
    return(summary_data[,1:4])
  }
  summary_data_sig <- summary_data[summary_data$p_val<alpha,]
  c <- stats::qnorm((alpha)/2, lower.tail=FALSE)
  beta.cl1 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl2 <- c(rep(0,nrow(summary_data_sig)))
  beta.cl3 <- c(rep(0,nrow(summary_data_sig)))
  
  for (i in 1:nrow(summary_data_sig)){
    cond.like <- function (mu) { return ((stats::dnorm(summary_data_sig$z[i]-mu))/(stats::pnorm(-c+mu)+stats::pnorm(-c-mu)))}
    mean.cond.like <- function (mu) {return(mu*cond.like(mu))}
    beta.cl1[i] <- (stats::optimize(cond.like, c(0,summary_data_sig$z[i]), maximum=TRUE)$maximum)*summary_data_sig$se[i]
    if(abs(summary_data_sig$z[i]) < 37){
      beta.cl2[i] <- ((stats::integrate(mean.cond.like,-37,37)$value)/(stats::integrate(cond.like,-37,37)$value))*(summary_data_sig$se[i])
    }else{
      if(abs(summary_data_sig$z[i]) > 100){
        beta.cl2[i] <- summary_data_sig$beta[i]
      }else{
        beta.cl2[i] <- ((stats::integrate(mean.cond.like,-100,100)$value)/(stats::integrate(cond.like,-100,100)$value))*(summary_data_sig$se[i])
      }
    }
    beta.cl3[i] <- (beta.cl1[i]+beta.cl2[i])/2
  }
  
  summary_data_sig <- cbind(summary_data_sig,beta.cl1,beta.cl2,beta.cl3)
  if(sum(abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) >= 1){
    bad.cl2 <- which((abs(beta.cl2) > abs(summary_data_sig$beta + sign(summary_data_sig$beta))) == TRUE)
    bad.cl2.df <- data.frame(rsid = summary_data_sig$rsid[bad.cl2], z = summary_data_sig$z[bad.cl2], z.cl2 = beta.cl2[bad.cl2]/summary_data_sig$se[bad.cl2], beta = summary_data_sig$beta[bad.cl2], beta.cl2 = beta.cl2[bad.cl2], se = summary_data_sig$se[bad.cl2])
    good.cl2.df <- data.frame(rsid = summary_data_sig$rsid[-(bad.cl2)], z = summary_data_sig$z[-(bad.cl2)], z.cl2 = beta.cl2[-(bad.cl2)]/summary_data_sig$se[-(bad.cl2)], beta = summary_data_sig$beta[-bad.cl2], beta.cl2 = beta.cl2[-(bad.cl2)], se = summary_data_sig$se[-(bad.cl2)])
    for(i in 1:nrow(bad.cl2.df)){
      if(length(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]]) == 0 || length(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]]) == 0){bad.cl2.df$z.cl2[i] <- good.cl2.df$z.cl2[which.min(abs(good.cl2.df$z-bad.cl2.df$z[i]))]
      bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]}else{
        min_greater <- min(good.cl2.df$z[good.cl2.df$z >= bad.cl2.df$z[i]])
        min_smaller <- min(good.cl2.df$z[good.cl2.df$z < bad.cl2.df$z[i]])
        bad.cl2.df$z.cl2[i] <- ((abs(min_greater - bad.cl2.df$z[i]))*(good.cl2.df$z.cl2[good.cl2.df$z == min_smaller]) + (abs(bad.cl2.df$z[i] -     min_smaller))*(good.cl2.df$z.cl2[good.cl2.df$z == min_greater]))/(abs(min_greater - min_smaller))
        bad.cl2.df$beta.cl2[i] <- bad.cl2.df$z.cl2[i]*bad.cl2.df$se[i]
      }
      summary_data_sig$beta.cl2[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- bad.cl2.df$beta.cl2[i]
      summary_data_sig$beta.cl3[summary_data_sig$rsid == bad.cl2.df$rsid[i]] <- (summary_data_sig$beta.cl1[summary_data_sig$rsid == bad.cl2.df$rsid[i]] +bad.cl2.df$beta.cl2[i])/2
    }
  }
  summary_data_sig <- dplyr::arrange(summary_data_sig,dplyr::desc(abs(summary_data_sig$z)))
  summary_data_sig <- cbind(summary_data_sig[,1:3],summary_data_sig[,6:8])
  return(summary_data_sig)
}

################################################################################
# 16. simulate_ss_bim() - Given values for heritability (H), polygenicity (Pi),
# sample size (nid) and selection coefficient (S), 'true' values of effect size
# and corresponding standard error are simulated in which effect sizes are
# assumed to have a bimodal distribution, i.e. 50% of effect sizes come from a
# normal distribution centered at 0 while the other half are generated from a
# normal distribution with mean 2.5.

simulate_ss_bim <- function(H,Pi,nid,sc,rep_nid=1){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- c(rnorm(0.5*effect_snps,2.5,sd=sqrt((2*maf*(1-maf))^sc)),rnorm(0.5*effect_snps,0,sd=sqrt((2*maf[(0.5*effect_snps+1):effect_snps]*(1-maf[(0.5*effect_snps+1):effect_snps]))^sc)))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  rep_se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid*nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se,rep_se)
  return(stats)
}

################################################################################
# 17. simulate_ss_exp() - Given values for heritability (H), polygenicity (Pi),
# sample size (nid) and selection coefficient (S), 'true' values of effect size
# and corresponding standard error are simulated in which effect sizes are
# assumed to have a skewed distribution.

simulate_ss_exp <- function(H,Pi,nid,sc,rep_nid=1){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- c(-rexp(n=.1*effect_snps,rate = 1/sqrt((2*maf*(1-maf))^sc)), rexp(n=.9*effect_snps,rate = 1/sqrt((2*maf[(0.1*effect_snps+1):effect_snps]*(1-maf[(0.1*effect_snps+1):effect_snps]))^sc)))
  var_y <- sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2)/H
  true_beta <- true_beta/sqrt(var_y)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(nid-2)*maf*(1-maf)))
  rep_se <- sqrt((1 - 2*maf*(1-maf)*true_beta^2)/(2*(rep_nid*nid-2)*maf*(1-maf)))
  stats <- data.frame(true_beta,se,rep_se)
  return(stats)
}

################################################################################
# 18. simulate_ss_bin() - Given values for heritability (H), polygenicity (Pi),
# sample size (nid) and selection coefficient (S), 'true' values of effect size
# and corresponding standard error are simulated for a binary trait in which
# effect sizes are assumed to have a normal distribution.

# Extra functions required first for logistic regression to obtain standard 
# errors:

# convert logit values to probabilities
expit <- function(x){exp(x)/(1+exp(x))}
# find Gamma0 for logistic model so prevalence is 0.1
myfunc <- function(MAF,Gamma0, Gamma1, prev=0.1){
  MAF^2*expit(Gamma0+Gamma1*2)+2*MAF*(1-MAF)*expit(Gamma0+Gamma1)+(1-MAF)^2*expit(Gamma0)-prev
}
# calculate asymptotic version of (X^t W X)
asymp_var_logistic <- function(n,G_prob,Gamma_0,Gamma_1){
  N <- length(Gamma_0)
  disease_probs <- expit(outer(Gamma_0,c(1,1,1)) +  outer(Gamma_1,c(0,1,2)))
  diag_weights <- disease_probs*(1-disease_probs)
  a <- n*apply(G_prob*diag_weights,1,sum)
  b <- n*apply(G_prob*diag_weights*matrix(rep(c(0,1,2),N),byrow=TRUE,nrow=N),1,sum)
  c <-  b
  d <- n*apply(G_prob*diag_weights*matrix(rep(c(0,1,4),N),byrow=TRUE,nrow=N),1,sum)
  # invert matrix and take element for Gamma1
  return(a/(a*d-b*c))
}
# grid-like computation of beta_0 corresponding to values for maf and true_beta
# with fixed prevalence of 0.1 saves computational time
params <- expand.grid(
  maf = seq(0.01,0.5,by=0.01),
  true_beta = seq(-1.5,1.5,by=0.01)
)
params$maf <- round(params$maf,2)
params$true_beta <- round(params$true_beta,2)
params$Gamma_0 <- c(rep(0,nrow(params)))
for(i in 1:nrow(params)) {params$Gamma_0[i] <- uniroot(myfunc,Gamma1=params$true_beta[i],MAF=params$maf[i],prev=0.1,lower=-10, upper=10)$root}

simulate_ss_bin <- function(H2,Pi,nid,sc){
  effect_snps <- Pi*n_snps
  maf <- runif(n_snps,0.01,0.5)
  true_beta <- rnorm(effect_snps,0,sd=sqrt((2*maf*(1-maf))^sc))
  scaling <- (1.6^2*H2)/((sum(2*maf[1:effect_snps]*(1-maf[1:effect_snps])*true_beta^2))*(1-H2))
  true_beta <- true_beta*sqrt(scaling)
  true_beta <- c(true_beta, rep(0,n_snps-effect_snps))
  rounded_maf <- round(maf,2)
  rounded_true_beta <- round(true_beta,2)
  Gamma0 <- c(rep(0,n_snps))
  dat <- data.frame(rounded_maf, rounded_true_beta)
  for(i in 1:n_snps) Gamma0[i] <- params$Gamma_0[7500 + 5000*rounded_true_beta[i] + 100*rounded_maf[i]]
  # assuming HWE
  G_prob <- cbind((1-maf)^2,2*maf*(1-maf),maf^2)
  Gamma1 <- true_beta
  var_Gamma_y <- asymp_var_logistic(nid,G_prob,Gamma0,Gamma1)
  se <- sqrt(var_Gamma_y)
  stats <- data.frame(true_beta,se)
  return(stats)
}

################################################################################
# 8. mse_sig_improve_per() - Given the data frame obtained due to method
# application, the true effect sizes and the significance threshold, the relative
# improvement in average MSE for significant SNPs due to method application is
# computed.

mse_sig_improve_per <- function(out,true_beta,i=5,alpha=5e-8){
  snp_sig <- out[abs(out$beta/out$se) > qnorm((alpha)/2, lower.tail=FALSE),]
  if (nrow(snp_sig) == 0){return(100)}
  mse_sig_improve <- (mean((true_beta[snp_sig$rsid]-snp_sig[,i])^2)-mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))/(mean((true_beta[snp_sig$rsid]-snp_sig$beta)^2))
  return(mse_sig_improve)
}
