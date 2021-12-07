# HMC simulations for the non-central 8 schools model
rm(list = ls())
set.seed(runif(1))

# library(dplyr)
library(doParallel)
registerDoParallel(cores = detectCores()-1)

# Importing half-t functions
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source('half_t/half_t_functions.R')
source('estimators.R')

# Metric considered: L2
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}

half_t_beta_l2 <- function(chain1, chain2){
  metric_l2(chain1$beta, chain2$beta)
}

# Import Riboflavin data
load("half_t/riboflavin.RData")
X <- as.matrix(riboflavin$x)
colnames(X) <- NULL
rownames(X) <- NULL
X_transpose <- t(X)
y <- as.matrix(riboflavin$y)
n <- dim(X)[1]
dimension <- dim(X)[2]

# # Generate synthetic data
# n <- 500
# dimension <- 50000
# s <- 20
# true_beta <- matrix(0,dimension,1)
# true_beta[1:s] = 2^(-(seq(s)/4-9/4))
# X <- matrix(rnorm(n*dimension), nrow = n, ncol = dimension)
# X_transpose <- t(X)
# #Error terms
# error_std <- 2
# error_terms = error_std*rnorm(n, mean = 0, sd = 1)
# y = X%*%true_beta + error_terms

# Simulation parameters
no_chains <- 100
burnin <- 1e3
chain_length <- 2e3 + burnin

approxMCMCparams <- c(0,0.0001,0.001,0.01)
bounds_df <- data.frame()
for (epsilon in approxMCMCparams){
  for(t_dist_df in c(2)){
    crn_half_t_sampler <- function(){
      coupled_chain <- 
        coupled_half_t_mcmc(mc_chain_size=chain_length, X, X_transpose, y, 
                            approximate_algo_delta_1=0, approximate_algo_delta_2=epsilon, 
                            a0=1, b0=1, std_MH=0.8, t_dist_df=t_dist_df, 
                            epsilon_eta=0, verbose=TRUE, metric=half_t_beta_l2, store_states=TRUE)
      return(list('P'=coupled_chain$beta_samples1[c((burnin+1):chain_length),],
                  'Q'=coupled_chain$beta_samples2[c((burnin+1):chain_length),]))
    }
    crn_bounds <- W2L2_UBLB_estimates(crn_half_t_sampler, no_chains=no_chains,
                                      lb='marginals', parallel=TRUE)
    bounds_df <-
      rbind(bounds_df,
            data.frame(n=n, dimension=dimension, approxMCMCparam=epsilon, t_dist_df=t_dist_df,
                       W2L2UBmean=crn_bounds$W2L2_UBmean, W2L2UBsd=crn_bounds$W2L2_UBsd,
                       W2L2LBmean=crn_bounds$W2L2_LBmean, W2L2LBsd=crn_bounds$W2L2_LBsd))
    
    print(c(t_dist_df, epsilon))
  }
}

# save(bounds_df, file="half_t/half_t_bounds_df_riboflavin.RData")
# save(bounds_df, file="half_t/half_t_bounds_df_synthetic.RData")
