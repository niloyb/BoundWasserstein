rm(list = ls())
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source(file = 'skinny_gibbs/skinny_gibbs_functions.R')
source(file = 'estimators.R')

library(foreach)
library(doParallel)
registerDoParallel(cores = detectCores()-1)

# Metrics considered: L2 and Hamming. 
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}
# Note: L2 distance squared now gives the Hamming distance
metric_hamming <- function(x,y){metric_l2(x,y)^2}

skinny_gibbs_df <- data.frame()

n <- 100
for(dimension in seq(100,500,100)){
  # Generate data
  s0 <- 4
  true_beta <- matrix(0,dimension,1)
  #true_beta[1:s0] = 2^(-(seq(s0)/4-9/4))
  true_beta[1:s0] = c(-1.5,2,-2.5,3) # true beta: following skinny gibbs paper
  X <- matrix(rnorm(n*dimension), nrow = n, ncol = dimension)
  # Error terms 
  true_aug_y = rlogis(n, location = X%*%true_beta) # recall variance is pi^2/3
  y <- ifelse(true_aug_y>0,1,0) # Logistic response
  X <- matrix(scale(X),n,dimension) ### important for Bayesian variable selection
  
  # Choice of q, tau0, tau1: following skinny gibbs paper
  K <- max(10,log(n))
  q_seq <- seq(0,1,0.0001)
  probs <- abs(pbinom(K,dimension,q_seq)-0.9)
  q_index <- which(probs==min(probs))
  q <- q_seq[q_index]
  tau0 <- 1/sqrt(n)
  tau1 <- 1 # sqrt(max(1, dimension^(2.1)/(100*n)))
  nu <- 7.3
  
  burnin <- 1e3
  chain_length <- burnin+1e3
  no_chains <- 10
  # # Exact chain
  # exact_chain <- full_gibbs(X, y, tau0, tau1, q, nu, chain_length=chain_length)
  # # matplot(exact_chain$beta[,c(1:10)], type = 'l')
  # # which(colMeans(exact_chain$z[(burnin:chain_length),])>0.5)
  # # Skinny chain
  # skinny_chain <- skinny_gibbs(X, y, tau0, tau1, q, nu, chain_length=chain_length)
  # # matplot(skinny_chain$beta[,c(1:10)], type = 'l')
  # # which(colMeans(skinny_chain$z[(burnin:chain_length),])>0.5)
  # # matplot(metric_hamming(exact_chain$z,skinny_chain$z), type = 'l')
  
  # CRN exact chain
  # crn_exact_chain <- coupled_exact_chain(X, y, tau0, tau1, q, nu, chain_length=chain_length)
  # matplot(metric_hamming(crn_exact_chain$z1,crn_exact_chain$z2), type = 'l')
  # CRN skinny chain
  # crn_skinny_chain <- coupled_skinny_chain(X, y, tau0, tau1, q, nu, chain_length=chain_length)
  # matplot(metric_hamming(crn_skinny_chain$z1,crn_skinny_chain$z2), type = 'l')
  
  # CRN coupling of skinny and exact z chain
  # crn_exact_skinny_chain <- exact_skinny_chain(X, y, tau0, tau1, q, nu=nu, chain_length=chain_length)
  # which(colMeans(crn_exact_skinny_chain$exact_z[(burnin:chain_length),])>0.5)
  # which(colMeans(crn_exact_skinny_chain$skinny_z[(burnin:chain_length),])>0.5)
  # matplot(metric_hamming(crn_exact_skinny_chain$exact_z,crn_exact_skinny_chain$skinny_z), type = 'l')
  
  crn_exact_skinny_chain_sampler_z <- function(){
    output <- exact_skinny_chain(X, y, tau0, tau1, q, nu=nu, chain_length=chain_length)
    return(list('P'=output$exact_z, 'Q'=output$skinny_z))
  }
  crn_bounds_z <- 
    W2L2_UBLB_estimates(crn_exact_skinny_chain_sampler_z, no_chains=no_chains, 
                        lb='marginals', parallel=FALSE)
  
  crn_exact_skinny_chain_sampler_beta <- function(){
    output <- exact_skinny_chain(X, y, tau0, tau1, q, nu=nu, chain_length=chain_length)
    return(list('P'=output$exact_beta, 'Q'=output$skinny_beta))
  }
  crn_bounds_beta <- 
    W2L2_UBLB_estimates(crn_exact_skinny_chain_sampler_beta, no_chains=no_chains, 
                        lb='marginals', parallel=FALSE)
  
  # Upper and lower bounds
  skinny_gibbs_df <- 
    rbind(skinny_gibbs_df,
          data.frame(n=n, dimension=dimension, sparsity=s0, no_chains=no_chains,
                     W2L2UBmean_z=crn_bounds_z$W2L2_UBmean, W2L2UBsd_z=crn_bounds_z$W2L2_UBsd,
                     W2L2LBmean_z=crn_bounds_z$W2L2_LBmean, W2L2LBsd_z=crn_bounds_z$W2L2_LBsd,
                     W2L2UBmean_beta=crn_bounds_beta$W2L2_UBmean, W2L2UBsd_beta=crn_bounds_beta$W2L2_UBsd,
                     W2L2LBmean_beta=crn_bounds_beta$W2L2_LBmean, W2L2LBsd_beta=crn_bounds_beta$W2L2_LBsd,
                     W2L2UBs_z=I(list(crn_bounds_z$W2L2_UBs)), W2L2LBs_z=I(list(crn_bounds_z$W2L2_LBs)),
                     W2L2UBs_beta=I(list(crn_bounds_beta$W2L2_UBs)), W2L2LBs_beta=I(list(crn_bounds_beta$W2L2_LBs))))
  # 
  # crn_exact_skinny_output <- 
  #   wp_ub_estimate(crn_exact_skinny_chain_sampler, no_chains = no_chains, p=2, metric = metric_l2)
  # 
  # coupled_traj <- crn_exact_skinny_output$wp_power_p_ub_mean_tracjectory
  # # metric_l2(crn_exact_skinny_chain$exact, crn_exact_skinny_chain$skinny)
  # coupled_W2L2UB <- mean(coupled_traj[burnin:chain_length])^(1/2)
  # coupled_W1hammingUBmean <- coupled_W2L2UB^2
  # coupled_W1hammingUBsd <- sd(rowMeans(crn_exact_skinny_output$wp_power_p_ub_tracjectories[,burnin:chain_length]))
  # 
  # indep_exact_skinny_chain_sampler <- function(){
  #   exact_chain <- full_gibbs(X, y, tau0, tau1, q, w, chain_length=chain_length)
  #   skinny_chain <- skinny_gibbs(X, y, tau0, tau1, q, w, chain_length=chain_length)
  #   return(list('P'=exact_chain$z, 'Q'=skinny_chain$z))
  # }
  # indep_exact_skinny_output <- 
  #   wp_ub_estimate(indep_exact_skinny_chain_sampler, no_chains = no_chains, p=2, metric = metric_l2)
  # indep_traj <- indep_exact_skinny_output$wp_power_p_ub_mean_tracjectory
  # # metric_l2(crn_exact_skinny_chain$exact, crn_exact_skinny_chain$skinny)
  # indep_W2L2UB <- mean(indep_traj[burnin:chain_length])^(1/2)
  # indep_W1hammingUBmean <- indep_W2L2UB^2
  # indep_W1hammingUBsd <- sd(rowMeans(indep_exact_skinny_output$wp_power_p_ub_tracjectories[,burnin:chain_length]))
  # 
  # W2L2_lbs <- w2l2_lb_estimate(indep_exact_skinny_chain_sampler, no_chains = no_chains)
  # W1hammingLBmean <- mean(W2L2_lbs$w2l2_lbs^2)
  # W1hammingLBsd <- sd(W2L2_lbs$w2l2_lbs^2)
  
  # # Upper and lower bounds
  # skinny_gibbs_df <- 
  #   rbind(skinny_gibbs_df,
  #         data.frame(n=n, p=p, sparsity=s0, no_chains=no_chains,
  #                    coupled_W1hammingUBmean=coupled_W1hammingUBmean, coupled_W1hammingUBsd, coupled_W1hammingUBsd,
  #                    W1hammingLBmean=W1hammingLBmean, W1hammingLBsd=W1hammingLBsd,
  #                    indep_W1hammingUBmean=indep_W1hammingUBmean, indep_W1hammingUBsd=indep_W1hammingUBsd))
  print(dimension)
}

# Saving data
# save(skinny_gibbs_df, file="skinny_gibbs/skinny_gibbs_bounds_df.RData")



## Misc code ##
# # Implementation from Skinny Gibbs package
# library(skinnybasad)
# B0=rep(0,p)
# ind = sample(seq(1: p), 5)
# Z0=rep(0,p); Z0[ind] = 1
# exact_chain2 = skinnybasad(X,E=y,pr=q,B0, Z0,nsplit=10,modif=0,nburn=burnin,niter=chain_length,printitrsep=20000,maxsize=30)
# which(exact_chain2$marZ>0.5) # Should be similar to own implementation
# skinny_chain2 = skinnybasad(X,E=y,pr=q,B0, Z0,nsplit=10,modif=0,nburn=burnin,niter=chain_length,printitrsep=20000,maxsize=30)
# which(skinny_chain2$marZ>0.5) # Should be similar to own implementation

# # CRN coupling of exact chain
# crn_exact_chain <- coupled_exact_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
# which(colMeans(crn_exact_chain$z1[(burnin:chain_length),])>0.5) # 
# which(colMeans(crn_exact_chain$z2[(burnin:chain_length),])>0.5) # 
# matplot((rowSums((crn_exact_chain$beta1-crn_exact_chain$beta2)^2)^0.5),
#         type = 'l', ylab = 'beta distance', log='y')
# matplot((rowSums((crn_exact_chain$z1-crn_exact_chain$z2)^2)^0.5),
#         type = 'l', ylab = 'z distance')
# 
# # CRN coupling of skinny chain
# crn_skinny_chain <- 
#   coupled_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
# which(colMeans(crn_skinny_chain$z1[(burnin:chain_length),])>0.5) # 
# which(colMeans(crn_skinny_chain$z2[(burnin:chain_length),])>0.5) # 
# matplot((rowSums((crn_skinny_chain$beta1-crn_skinny_chain$beta2)^2)^0.5),
#         type = 'l', ylab = 'beta distance', log='y')
# matplot((rowSums((crn_skinny_chain$z1-crn_skinny_chain$z2)^2)^0.5),
#         type = 'l', ylab = 'z distance')


