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

######################## UCI Malware dataset ########################
malware_data <- read.csv('/Users/niloybiswas/Downloads/uci_malware_detection.csv', header = TRUE)
y <- as.matrix(rep(0, nrow(malware_data)))
y[malware_data[,1]=='malicious',] <- 1
X <- as.matrix(malware_data[,-1])
rownames(X) <- NULL
colnames(X) <- NULL

# Removing covariates which do not vary across observations
col_sds <- sapply(c(1:ncol(X)), function(i){sd(X[,i])})
X <- X[,col_sds!=0]
n <- nrow(X)
dimension <- ncol(X)
X <- as.matrix(scale(X),n,dimension)

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
no_chains <- 100

######## W2L2 upper and lower bounds ########
crn_exact_skinny_chain_sampler_beta <- function(){
  output <- exact_skinny_chain(X, y, tau0, tau1, q, nu=nu, chain_length=chain_length)
  return(list('P'=output$exact_beta[-c(1:burnin),], 'Q'=output$skinny_beta[-c(1:burnin),]))
}
crn_bounds_beta <- 
  W2L2_UBLB_estimates(crn_exact_skinny_chain_sampler_beta, no_chains=no_chains, 
                      lb='marginals', parallel=FALSE)
# Upper and lower bounds
skinny_gibbs_df_malware <- 
  data.frame(dataset='malware', n=n, dimension=dimension, no_chains=no_chains,
             W2L2UBmean_beta=crn_bounds_beta$W2L2_UBmean, W2L2UBsd_beta=crn_bounds_beta$W2L2_UBsd,
             W2L2LBmean_beta=crn_bounds_beta$W2L2_LBmean, W2L2LBsd_beta=crn_bounds_beta$W2L2_LBsd,
             W2L2UBs_beta=I(list(crn_bounds_beta$W2L2_UBs)), W2L2LBs_beta=I(list(crn_bounds_beta$W2L2_LBs)))


######################## Lymph Node dataset ########################
lymph_data <- read.table('skinny_gibbs/lymph.dat', header = FALSE)
# dim(lymph_data)

y <- lymph_data[,1]
X <- lymph_data[,-1]
colnames(X) <- NULL
rownames(X) <- NULL 
n <- nrow(X)
dimension <- ncol(X)
X <- as.matrix(scale(X),n,dimension)

# Choice of q, tau0, tau1: following skinny Gibbs paper
K <- max(10,log(n))
q_seq <- seq(0,1,0.0001)
probs <- abs(pbinom(K,dimension,q_seq)-0.9)
q_index <- which(probs==min(probs))
q <- q_seq[q_index]
tau0 <- 1/sqrt(n)
tau1 <- 1 # sqrt(max(1, dimension^(2.1)/(100*n)))

# Fixing w
s <- pi/sqrt(3) # sd of logistic distribution
w <- rep(1,n)/(s^2)

burnin <- 1e3
chain_length <- burnin+2e3
no_chains <- 100

######## W2L2 upper and lower bounds ########
crn_exact_skinny_chain_sampler_beta <- function(){
  output <- exact_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
  return(list('P'=output$exact_beta[-c(1:burnin),], 'Q'=output$skinny_beta[-c(1:burnin),]))
}
crn_bounds_beta <- 
  W2L2_UBLB_estimates(crn_exact_skinny_chain_sampler_beta, no_chains=no_chains, 
                      lb='marginals', parallel=FALSE)
# Upper and lower bounds
skinny_gibbs_df_lymph <- 
  data.frame(dataset='lymph', n=n, dimension=dimension, no_chains=no_chains,
             W2L2UBmean_beta=crn_bounds_beta$W2L2_UBmean, W2L2UBsd_beta=crn_bounds_beta$W2L2_UBsd,
             W2L2LBmean_beta=crn_bounds_beta$W2L2_LBmean, W2L2LBsd_beta=crn_bounds_beta$W2L2_LBsd,
             W2L2UBs_beta=I(list(crn_bounds_beta$W2L2_UBs)), W2L2LBs_beta=I(list(crn_bounds_beta$W2L2_LBs)))


skinny_gibbs_df <- rbind(skinny_gibbs_df_malware, skinny_gibbs_df_lymph)
# Saving data
# save(skinny_gibbs_df, file="skinny_gibbs/skinny_gibbs_df.RData")

