# Logistic regression plots

rm(list = ls())

# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(rstan)

library(doParallel)
registerDoParallel(cores = detectCores()-1)

# Functions
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source('kernels.R')
source('estimators.R')
source('rcpp_functions.R')
source('bayesian_logistic_regression/logreg_functions.R')

# Metric considered: L2
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}

# Wasserstein's distance between Gaussians
W2L2mvns <- function(mu1,mu2,Sigma1,Sigma2){
  sqrtSigma1 <- expm::sqrtm(Sigma1)
  crossSigma <- expm::sqrtm(sqrtSigma1%*%(Sigma2%*%sqrtSigma1))
  out <- sum(diag(Sigma1+Sigma2-2*crossSigma)) + sum((mu1-mu2)^2)
  return(out^0.5)
}

# p Wasserstein's distance
p <- 2

# Dataset (synthetic, pima, ds1 considered)
dataset <- 'ds1'
if (dataset == 'synthetic'){
  n <- 1000
  dimension <- 10
  X <- matrix(rnorm(n*dimension), nrow = n, ncol = dimension)
  s <- 5
  true_beta <- matrix(0,dimension,1)
  true_beta[1:s] = 2^(-(seq(s)/4-9/4))
  probs <- exp(logSigmoid(X%*%true_beta))
  y <- 2*rbinom(n, size = 1, prob = probs)-1
} else if (dataset == 'pima'){
  pima_indians_csv <-
    read.csv('bayesian_logistic_regression/pima-indians-diabetes-new.csv', header = TRUE)
  y <- (2*pima_indians_csv[,9]-1)
  X <- as.matrix(pima_indians_csv[,c(1:8)])
} else if (dataset == 'ds1'){
  ds1_csv <- read.csv("bayesian_logistic_regression/ds1.10.csv.bz2")
  y <- 2*ds1_csv[,11]-1
  X <- as.matrix(ds1_csv[,c(1:10)])
}
colnames(X) <- NULL
X <- scale(X) # standardising the design matrix
n <- dim(X)[1]
dimension <- dim(X)[2]

# Prior
prior_mean <- rep(0,dimension)
inv_prior_cov_matrix <- diag(dimension)

############################# Simulating data ##################################
if ((dataset == 'synthetic') | (dataset == 'pima')){
  # longer chain and smaller step size for ds1 dataset where n<=1000
  burnin <- 1e3
  chain_length <- 1e3+burnin
  no_chains <- 100
  sigma_mh <- 0.05
} else if (dataset == 'ds1'){
  # shorter chain and smaller step size for ds1 dataset where n>=20000
  burnin <- 500
  chain_length <- 500+burnin
  no_chains <- 40
  sigma_mh <- 0.01 
}

############################# ADVI (mean field) ################################
# Compiling Stan model
data_stan <- 
  list('n'=n, 'dimension' = dimension, 'x' = X, 'y'=(y+1)/2, 'prior_mean'=prior_mean, 
       'prior_cov'=solve(inv_prior_cov_matrix))
logreg_stan <- stan_model(file = "bayesian_logistic_regression/logreg_model.stan")
mf_ADVI <- vb(logreg_stan, data=data_stan, output_samples=10000, seed=1,
              tol_rel_obj=0.0001)
advi_output <- extract(mf_ADVI)
mu <- colMeans(advi_output$beta)
inv_sigma <- diag(1/apply(advi_output$beta,2,var))
mu_ADVI <- mu
sigma_ADVI <- solve(inv_sigma)
coupled_chain_ADVI <- function(){
  # Starting beta1_init from (approximately) stationarity
  # NOTE: can use L-Lag Couplings for principles burn-in guidance
  beta1_init <- rnorm(dimension)
  beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
  beta1_init <- beta1_long_run[2e3,]
  beta2_init <- mvtnorm::rmvnorm(1, mean=mu_ADVI, sigma=sigma_ADVI)[1,]
  # beta2_init <- rnorm(dimension)
  coupled_mala(beta1_init, beta2_init, targetLogPdf, 
               mvnLogPdf, targetGradLogPdf, mvnGradLogPdf, sigma_mh, chain_length,
               reflect_threshold=0)
}

crn_mf_ADVI <- W2L2_UBLB_estimates(coupled_chain_ADVI, no_chains=no_chains)
crn_mf_W2L2UBmean <- crn_mf_ADVI$W2L2_UBmean
crn_mf_W2L2UBsd <- crn_mf_ADVI$W2L2_UBsd
mf_W2L2LBmean <- crn_mf_ADVI$W2L2_LBmean
mf_W2L2LBsd <- crn_mf_ADVI$W2L2_LBsd

# crn_mf_ADVI <- wp_ub_estimate(coupled_chain_ADVI, no_chains=no_chains,
#                               p=p, metric=metric_l2, parallel=TRUE)
# crn_mf_traj <- crn_mf_ADVI$wp_power_p_ub_mean_tracjectory
# crn_mf_W2L2UBmean <- mean(crn_mf_traj[burnin:chain_length])^(1/2)
# crn_mf_W2L2UBsd <- sd(rowMeans(crn_mf_ADVI$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))

# indep_mf_ADVI_sampler <- function(){
#   # Starting beta1_init from (approximately) stationarity
#   # NOTE: can use L-Lag Couplings for principles burn-in guidance
#   beta1_init <- rnorm(dimension)
#   beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
#   beta1_init <- beta1_long_run[2e3,]
#   beta1_chain <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, chain_length)
#   beta2_init <- mvtnorm::rmvnorm(1, mean=mu_ADVI, sigma=sigma_ADVI)[1,]
#   beta2_chain <- mala(beta2_init, mvnLogPdf, mvnGradLogPdf, sigma_mh, chain_length)
#   return(list('P'=beta1_chain, 'Q'=beta2_chain))
# }
# indep_mf_ADVI_output <- 
#   wp_ub_estimate(indep_mf_ADVI_sampler, no_chains = no_chains, p=2, metric = metric_l2)
# indep_traj <- indep_mf_ADVI_output$wp_power_p_ub_mean_tracjectory
# # metric_l2(crn_exact_skinny_chain$exact, crn_exact_skinny_chain$skinny)
# indep_mf_W2L2UBmean <- mean(indep_traj[burnin:chain_length])^(1/2)
# indep_mf_W2L2UBsd <- sd(rowMeans(indep_mf_ADVI_output$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))

# mf_W2L2_lbs <- w2l2_lb_estimate(indep_mf_ADVI_sampler, no_chains = no_chains)
# mf_W2L2LBmean <- mean(mf_W2L2_lbs$w2l2_lbs)
# mf_W2L2LBsd <- sd(mf_W2L2_lbs$w2l2_lbs)

################ Laplace Approximation ################
logreg_laplace <- optim(rep(0,dimension), method='BFGS', targetLogPdf, control = list(fnscale=-1), hessian = TRUE)
mu <- logreg_laplace$par
inv_sigma <- -logreg_laplace$hessian
mu_laplace <- mu
sigma_laplace <- solve(inv_sigma)
coupled_chain_laplace <- function(){
  # Starting beta1_init from (approximately) stationarity
  # NOTE: can use L-Lag Couplings to do this in a more principled way
  beta1_init <- rnorm(dimension)
  beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
  beta1_init <- beta1_long_run[2e3,]
  beta2_init <- mvtnorm::rmvnorm(1, mean=mu_laplace, sigma=sigma_laplace)[1,]
  # beta2_init <- rnorm(dimension)
  coupled_mala(beta1_init, beta2_init, targetLogPdf, 
               mvnLogPdf, targetGradLogPdf, mvnGradLogPdf, sigma_mh, chain_length,
               reflect_threshold=0)
}

crn_laplace <- W2L2_UBLB_estimates(coupled_chain_laplace, no_chains=no_chains)
crn_laplace_W2L2UBmean <- crn_laplace$W2L2_UBmean
crn_laplace_W2L2UBsd <- crn_laplace$W2L2_UBsd
laplace_W2L2LBmean <- crn_laplace$W2L2_LBmean
laplace_W2L2LBsd <- crn_laplace$W2L2_LBsd





# crn_laplace <- wp_ub_estimate(coupled_chain_laplace, no_chains=no_chains,
#                               p=p, metric=metric_l2, parallel=TRUE)
# crn_laplace_traj <- crn_laplace$wp_power_p_ub_mean_tracjectory
# crn_laplace_W2L2UBmean <- mean(crn_laplace_traj[burnin:chain_length])^(1/2)
# crn_laplace_W2L2UBsd <- sd(rowMeans(crn_laplace$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# indep_laplace_sampler <- function(){
#   # Starting beta1_init from (approximately) stationarity
#   # NOTE: can use L-Lag Couplings for principles burn-in guidance
#   beta1_init <- rnorm(dimension)
#   beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
#   beta1_init <- beta1_long_run[2e3,]
#   beta1_chain <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, chain_length)
#   beta2_init <- mvtnorm::rmvnorm(1, mean=mu_ADVI, sigma=sigma_ADVI)[1,]
#   beta2_chain <- mala(beta2_init, mvnLogPdf, mvnGradLogPdf, sigma_mh, chain_length)
#   return(list('P'=beta1_chain, 'Q'=beta2_chain))
# }
# indep_laplace_output <- 
#   wp_ub_estimate(indep_laplace_sampler, no_chains = no_chains, p=2, metric = metric_l2)
# indep_traj <- indep_laplace_output$wp_power_p_ub_mean_tracjectory
# indep_laplace_W2L2UBmean <- mean(indep_traj[burnin:chain_length])^(1/2)
# indep_laplace_W2L2UBsd <- sd(rowMeans(indep_laplace_output$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# laplace_W2L2_lbs <- w2l2_lb_estimate(indep_laplace_sampler, no_chains = no_chains)
# laplace_W2L2LBmean <- mean(laplace_W2L2_lbs$w2l2_lbs)
# laplace_W2L2LBsd <- sd(laplace_W2L2_lbs$w2l2_lbs)

################ SGLD 0.1 ################
subsample_ratio <- 0.1
coupled_chain_sgld0.1 <- function(){
  # Starting beta1_init from (approximately) stationarity
  # NOTE: can use L-Lag Couplings for more principled burn-in guidance
  beta1_init <- rnorm(dimension)
  beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
  beta1_init <- beta1_long_run[2e3,]
  beta2_init <- rnorm(dimension)
  beta2_long_run <- ula(beta2_init, approxTargetGradLogPdf, sigma_mh, 2e3)
  beta2_init <- beta2_long_run[2e3,]
  coupled_mala(beta1_init, beta2_init, targetLogPdf, NA, 
               targetGradLogPdf, approxTargetGradLogPdf, sigma_mh, chain_length,
               reflect_threshold=0, pi_correction = TRUE, q_correction = FALSE)
}
sgld0.1_chain <- ula(rnorm(dimension), approxTargetGradLogPdf, sigma_mh, 1e4)
mu_sgld0.1 <- colMeans(sgld0.1_chain[c(5e3:1e4),])
sigma_sgld0.1 <- cov(sgld0.1_chain[c(5e3:1e4),])

crn_sgld0.1 <- W2L2_UBLB_estimates(coupled_chain_sgld0.1, no_chains=no_chains)
crn_sgld0.1_W2L2UBmean <- crn_sgld0.1$W2L2_UBmean
crn_sgld0.1_W2L2UBsd <- crn_sgld0.1$W2L2_UBsd
sgld0.1_W2L2LBmean <- crn_sgld0.1$W2L2_LBmean
sgld0.1_W2L2LBsd <- crn_sgld0.1$W2L2_LBsd

# crn_sgld0.1 <- wp_ub_estimate(coupled_chain_sgld0.1, no_chains=no_chains, p=p,
#                               metric=metric_l2, parallel=TRUE)
# crn_sgld0.1_traj <- crn_sgld0.1$wp_power_p_ub_mean_tracjectory
# crn_sgld0.1_W2L2UBmean <- mean(crn_sgld0.1_traj[burnin:chain_length])^(1/2)
# crn_sgld0.1_W2L2UBsd <- sd(rowMeans(crn_sgld0.1$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# indep_sgld0.1_sampler <- function(){
#   # Starting beta1_init from (approximately) stationarity
#   # NOTE: can use L-Lag Couplings for principles burn-in guidance
#   beta1_init <- rnorm(dimension)
#   beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
#   beta1_init <- beta1_long_run[2e3,]
#   beta1_chain <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, chain_length)
#   beta2_init <- rnorm(dimension)
#   beta2_long_run <- ula(beta2_init, approxTargetGradLogPdf, sigma_mh, 2e3)
#   beta2_init <- beta2_long_run[2e3,]
#   beta2_chain <- ula(beta2_init, approxTargetGradLogPdf, sigma_mh, chain_length)
#   return(list('P'=beta1_chain, 'Q'=beta2_chain))
# }
# indep_sgld0.1_output <- 
#   wp_ub_estimate(indep_sgld0.1_sampler, no_chains = no_chains, p=2, metric = metric_l2)
# indep_traj <- indep_sgld0.1_output$wp_power_p_ub_mean_tracjectory
# indep_sgld0.1_W2L2UBmean <- mean(indep_traj[burnin:chain_length])^(1/2)
# indep_sgld0.1_W2L2UBsd <- sd(rowMeans(indep_sgld0.1_output$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# sgld0.1_W2L2_lbs <- w2l2_lb_estimate(indep_sgld0.1_sampler, no_chains = no_chains)
# sgld0.1_W2L2LBmean <- mean(sgld0.1_W2L2_lbs$w2l2_lbs)
# sgld0.1_W2L2LBsd <- sd(sgld0.1_W2L2_lbs$w2l2_lbs)

################ SGLD 0.5 ################
subsample_ratio <- 0.5
coupled_chain_sgld0.5 <- function(){
  # Starting beta1_init from (approximately) stationarity
  # NOTE: can use L-Lag Couplings for principles burn-in guidance
  beta1_init <- rnorm(dimension)
  beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
  beta1_init <- beta1_long_run[2e3,]
  beta2_init <- rnorm(dimension)
  beta2_long_run <- ula(beta2_init, approxTargetGradLogPdf, sigma_mh, 2e3)
  beta2_init <- beta2_long_run[2e3,]
  coupled_mala(beta1_init, beta2_init, targetLogPdf, NA, 
               targetGradLogPdf, approxTargetGradLogPdf, sigma_mh, chain_length,
               reflect_threshold=0, pi_correction = TRUE, q_correction = FALSE)
}
sgld0.5_chain <- ula(rnorm(dimension), approxTargetGradLogPdf, sigma_mh, 1e4)
mu_sgld0.5 <- colMeans(sgld0.5_chain[c(5e3:1e4),])
sigma_sgld0.5 <- cov(sgld0.5_chain[c(5e3:1e4),])

crn_sgld0.5 <- W2L2_UBLB_estimates(coupled_chain_sgld0.5, no_chains=no_chains)
crn_sgld0.5_W2L2UBmean <- crn_sgld0.5$W2L2_UBmean
crn_sgld0.5_W2L2UBsd <- crn_sgld0.5$W2L2_UBsd
sgld0.5_W2L2LBmean <- crn_sgld0.5$W2L2_LBmean
sgld0.5_W2L2LBsd <- crn_sgld0.5$W2L2_LBsd


# crn_sgld0.5 <- wp_ub_estimate(coupled_chain_sgld0.5, no_chains=no_chains, p=p,
#                               metric=metric_l2, parallel=TRUE)
# crn_sgld0.5_traj <- crn_sgld0.5$wp_power_p_ub_mean_tracjectory
# crn_sgld0.5_W2L2UBmean <- mean(crn_sgld0.5_traj[burnin:chain_length])^(1/2)
# crn_sgld0.5_W2L2UBsd <- 
#   sd(rowMeans(crn_sgld0.5$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# indep_sgld0.5_sampler <- function(){
#   # Starting beta1_init from (approximately) stationarity
#   # NOTE: can use L-Lag Couplings for principles burn-in guidance
#   beta1_init <- rnorm(dimension)
#   beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
#   beta1_init <- beta1_long_run[2e3,]
#   beta1_chain <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, chain_length)
#   beta2_init <- rnorm(dimension)
#   beta2_long_run <- ula(beta2_init, approxTargetGradLogPdf, sigma_mh, 2e3)
#   beta2_init <- beta2_long_run[2e3,]
#   beta2_chain <- ula(beta2_init, approxTargetGradLogPdf, sigma_mh, chain_length)
#   return(list('P'=beta1_chain, 'Q'=beta2_chain))
# }
# indep_sgld0.5_output <- 
#   wp_ub_estimate(indep_sgld0.5_sampler, no_chains = no_chains, p=2, metric = metric_l2)
# indep_traj <- indep_sgld0.5_output$wp_power_p_ub_mean_tracjectory
# indep_sgld0.5_W2L2UBmean <- mean(indep_traj[burnin:chain_length])^(1/2)
# indep_sgld0.5_W2L2UBsd <- sd(rowMeans(indep_sgld0.5_output$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# sgld0.5_W2L2_lbs <- w2l2_lb_estimate(indep_sgld0.5_sampler, no_chains = no_chains)
# sgld0.5_W2L2LBmean <- mean(sgld0.5_W2L2_lbs$w2l2_lbs)
# sgld0.5_W2L2LBsd <- sd(sgld0.5_W2L2_lbs$w2l2_lbs)



################ ULA ################
coupled_chain_ula <- function(){
  # Starting beta1_init from (approximately) stationarity
  # NOTE: can use L-Lag Couplings for principles burn-in guidance
  beta1_init <- rnorm(dimension)
  beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
  beta1_init <- beta1_long_run[2e3,]
  beta2_init <- rnorm(dimension)
  beta2_long_run <- ula(beta2_init, targetGradLogPdf, sigma_mh, 2e3)
  beta2_init <- beta2_long_run[2e3,]
  coupled_mala(beta1_init, beta2_init, targetLogPdf, 
               NA, targetGradLogPdf, targetGradLogPdf, sigma_mh, chain_length,
               reflect_threshold=0, pi_correction = TRUE, q_correction = FALSE)
}
ula_chain <- ula(rnorm(dimension), targetGradLogPdf, sigma_mh, 1e4)
mu_ula <- colMeans(ula_chain[c(5e3:1e4),])
sigma_ula <- cov(ula_chain[c(5e3:1e4),])

crn_ula <- W2L2_UBLB_estimates(coupled_chain_ula, no_chains=no_chains)
crn_ula_W2L2UBmean <- crn_ula$W2L2_UBmean
crn_ula_W2L2UBsd <- crn_ula$W2L2_UBsd
ula_W2L2LBmean <- crn_ula$W2L2_LBmean
ula_W2L2LBsd <- crn_ula$W2L2_LBsd


# crn_ula <- wp_ub_estimate(coupled_chain_ula, no_chains=no_chains, 
#                           p=p, metric=metric_l2, parallel=TRUE)
# crn_ula_traj <- crn_ula$wp_power_p_ub_mean_tracjectory
# crn_ula_W2L2UBmean <- mean(crn_ula_traj[burnin:chain_length])^(1/2)
# crn_ula_W2L2UBsd <- sd(rowMeans(crn_ula$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# indep_ula_sampler <- function(){
#   # Starting beta1_init from (approximately) stationarity
#   # NOTE: can use L-Lag Couplings for principles burn-in guidance
#   beta1_init <- rnorm(dimension)
#   beta1_long_run <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 2e3)
#   beta1_init <- beta1_long_run[2e3,]
#   beta1_chain <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, chain_length)
#   beta2_init <- rnorm(dimension)
#   beta2_long_run <- ula(beta2_init, targetGradLogPdf, sigma_mh, 2e3)
#   beta2_init <- beta2_long_run[2e3,]
#   beta2_chain <- ula(beta2_init, targetGradLogPdf, sigma_mh, chain_length)
#   return(list('P'=beta1_chain, 'Q'=beta2_chain))
# }
# indep_ula_output <- 
#   wp_ub_estimate(indep_ula_sampler, no_chains = no_chains, p=2, metric = metric_l2)
# indep_traj <- indep_ula_output$wp_power_p_ub_mean_tracjectory
# indep_ula_W2L2UBmean <- mean(indep_traj[burnin:chain_length])^(1/2)
# indep_ula_W2L2UBsd <- sd(rowMeans(indep_ula_output$wp_power_p_ub_tracjectories[,burnin:chain_length])^(1/2))
# 
# ula_W2L2_lbs <- w2l2_lb_estimate(indep_ula_sampler, no_chains = no_chains)
# ula_W2L2LBmean <- mean(ula_W2L2_lbs$w2l2_lbs)
# ula_W2L2LBsd <- sd(ula_W2L2_lbs$w2l2_lbs)


################ Saving simulated data NEW ################
bounds_df <- 
  rbind(data.frame(dataset=dataset, type='mf_advi', n=n, dimension=dimension, stepsize=sigma_mh, 
                   no_chains=no_chains, W2L2UBmean=crn_mf_W2L2UBmean, W2L2UBsd=crn_mf_W2L2UBsd,
                   W2L2LBmean=mf_W2L2LBmean, W2L2LBsd=mf_W2L2LBsd),
        data.frame(dataset=dataset, type='laplace', n=n, dimension=dimension, stepsize=sigma_mh, 
                   no_chains=no_chains, W2L2UBmean=crn_laplace_W2L2UBmean, W2L2UBsd=crn_laplace_W2L2UBsd,
                   W2L2LBmean=laplace_W2L2LBmean, W2L2LBsd=laplace_W2L2LBsd),
        data.frame(dataset=dataset, type='sgld0.1', n=n, dimension=dimension, stepsize=sigma_mh, 
                   no_chains=no_chains, W2L2UBmean=crn_sgld0.1_W2L2UBmean, W2L2UBsd=crn_sgld0.1_W2L2UBsd,
                   W2L2LBmean=sgld0.1_W2L2LBmean, W2L2LBsd=sgld0.1_W2L2LBsd),
        data.frame(dataset=dataset, type='sgld0.5', n=n, dimension=dimension, stepsize=sigma_mh, 
                   no_chains=no_chains, W2L2UBmean=crn_sgld0.5_W2L2UBmean, W2L2UBsd=crn_sgld0.5_W2L2UBsd,
                   W2L2LBmean=sgld0.5_W2L2LBmean, W2L2LBsd=sgld0.5_W2L2LBsd),
        data.frame(dataset=dataset, type='ula', n=n, dimension=dimension, stepsize=sigma_mh, 
                   no_chains=no_chains, W2L2UBmean=crn_ula_W2L2UBmean, W2L2UBsd=crn_ula_W2L2UBsd,
                   W2L2LBmean=ula_W2L2LBmean, W2L2LBsd=ula_W2L2LBsd))

# filename <- paste("bayesian_logistic_regression/logreg_bounds_df2_",dataset,".RData",sep="")
# save(bounds_df, file=filename)









# ################ Saving simulated data ################
# # Trajectories
# trajectory <-
#   rbind(data.frame(t=c(1:chain_length), dataset=dataset, n=n, dimension=dimension, type='mf_advi',
#                    stepsize=sigma_mh, metric=(crn_mf_ADVI$wp_power_p_ub_mean_tracjectory)),
#         data.frame(t=c(1:chain_length), dataset=dataset, n=n, dimension=dimension, type='laplace', 
#                    stepsize=sigma_mh, metric=(crn_laplace$wp_power_p_ub_mean_tracjectory)),
#         data.frame(t=c(1:chain_length), dataset=dataset,  n=n, dimension=dimension, type='sgld0.1',
#                    stepsize=sigma_mh, metric=(crn_sgld0.1$wp_power_p_ub_mean_tracjectory)),
#         data.frame(t=c(1:chain_length), dataset=dataset, n=n, dimension=dimension, type='sgld0.5',
#                    stepsize=sigma_mh, metric=(crn_sgld0.5$wp_power_p_ub_mean_tracjectory)),
#         data.frame(t=c(1:chain_length), dataset=dataset, n=n, dimension=dimension, type='ula', 
#                    stepsize=sigma_mh, metric=(crn_ula$wp_power_p_ub_mean_tracjectory)))
# 
# # Marginal errors
# mala_long_run <- mala(rnorm(dimension), targetLogPdf, targetGradLogPdf, sigma_mh, 1e4)
# mu_mala <- colMeans(mala_long_run[c(5e3:1e4),])
# sigma_mala <- cov(mala_long_run[c(5e3:1e4),])
# 
# marginal_errors <-
#   data.frame('approx'=c("mf_advi", "laplace", "sgld0.1", "sgld0.5", "ula"),
#              dataset=dataset,
#              'relaive_cov_error'=
#                c(norm(solve(sigma_mala)%*%sigma_ADVI-diag(dimension), type = '2'),
#                  norm(solve(sigma_mala)%*%sigma_laplace-diag(dimension), type = '2'),
#                  norm(solve(sigma_mala)%*%sigma_sgld0.1-diag(dimension), type = '2'),
#                  norm(solve(sigma_mala)%*%sigma_sgld0.5-diag(dimension), type = '2'),
#                  norm(solve(sigma_mala)%*%sigma_ula-diag(dimension), type = '2')),
#              'cov_error'=
#                c(norm(solve(sigma_mala)-sigma_ADVI, type = '2'),
#                  norm(solve(sigma_mala)-sigma_laplace, type = '2'),
#                  norm(solve(sigma_mala)-sigma_sgld0.1, type = '2'),
#                  norm(solve(sigma_mala)-sigma_sgld0.5, type = '2'),
#                  norm(solve(sigma_mala)-sigma_ula, type = '2')),
#              'W2L2UB'=
#                c(crn_mf_ADVI$wp_ub,crn_laplace$wp_ub, crn_sgld0.1$wp_ub,
#                  crn_sgld0.5$wp_ub, crn_ula$wp_ub),
#              'W2L2LB'=
#                c(W2L2mvns(mu_mala, mu_ADVI, sigma_mala, sigma_ADVI),
#                  W2L2mvns(mu_mala, mu_laplace, sigma_mala, sigma_laplace),
#                  W2L2mvns(mu_mala, mu_sgld0.1, sigma_mala, sigma_sgld0.1),
#                  W2L2mvns(mu_mala, mu_sgld0.5, sigma_mala, sigma_sgld0.5),
#                  W2L2mvns(mu_mala, mu_ula, sigma_mala, sigma_ula)))
# 
# logreg_df <- list(trajectory = trajectory, marginal_errors = marginal_errors)
# 
# # filename <- paste("bayesian_logistic_regression/logreg_trajectory_df_",dataset,".RData",sep="")
# # save(logreg_df, file=filename)





