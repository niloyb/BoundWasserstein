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


######################## Malware dataset ########################
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
no_chains <- 10

######## W2L2 upper and lower bounds ########
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
skinny_gibbs_genome_df_malware <- 
  data.frame(dataset='malware', n=n, dimension=dimension, no_chains=no_chains,
             W2L2UBmean_z=crn_bounds_z$W2L2_UBmean, W2L2UBsd_z=crn_bounds_z$W2L2_UBsd,
             W2L2LBmean_z=crn_bounds_z$W2L2_LBmean, W2L2LBsd_z=crn_bounds_z$W2L2_LBsd,
             W2L2UBmean_beta=crn_bounds_beta$W2L2_UBmean, W2L2UBsd_beta=crn_bounds_beta$W2L2_UBsd,
             W2L2LBmean_beta=crn_bounds_beta$W2L2_LBmean, W2L2LBsd_beta=crn_bounds_beta$W2L2_LBsd,
             W2L2UBs_z=I(list(crn_bounds_z$W2L2_UBs)), W2L2LBs_z=I(list(crn_bounds_z$W2L2_LBs)),
             W2L2UBs_beta=I(list(crn_bounds_beta$W2L2_UBs)), W2L2LBs_beta=I(list(crn_bounds_beta$W2L2_LBs)))

######################## Lymph dataset ########################
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
chain_length <- burnin+1e3
no_chains <- 10

# # Exact chain
# exact_chain <- full_gibbs(X, y, tau0, tau1, q, w, chain_length=chain_length)
# which(colMeans(exact_chain$z[(burnin:chain_length),])>0.5)
# matplot(exact_chain$beta[(burnin:chain_length),colMeans(exact_chain$z[(burnin:chain_length),])>0.5], type='l')
# # Skinny chain
# skinny_chain <- skinny_gibbs(X, y, tau0, tau1, q, w, chain_length=chain_length)
# # which(colMeans(skinny_chain$z[(burnin:chain_length),])>0.5)
# matplot(skinny_chain$beta[(burnin:chain_length),colMeans(skinny_chain$z[(burnin:chain_length),])>0.5], type='l')

# crn_exact_chain <- coupled_exact_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
# crn_skinny_chain <- coupled_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
# matplot(metric_hamming(crn_exact_chain$z1,crn_exact_chain$z2), type = 'l')
# matplot(metric_hamming(crn_skinny_chain$z1,crn_skinny_chain$z2), type = 'l')


######## W2L2 upper and lower bounds ########
crn_exact_skinny_chain_sampler_z <- function(){
  output <- exact_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
  return(list('P'=output$exact_z, 'Q'=output$skinny_z))
}
crn_bounds_z <- 
  W2L2_UBLB_estimates(crn_exact_skinny_chain_sampler_z, no_chains=no_chains, 
                      lb='marginals', parallel=FALSE)

crn_exact_skinny_chain_sampler_beta <- function(){
  output <- exact_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
  return(list('P'=output$exact_beta, 'Q'=output$skinny_beta))
}
crn_bounds_beta <- 
  W2L2_UBLB_estimates(crn_exact_skinny_chain_sampler_beta, no_chains=no_chains, 
                      lb='marginals', parallel=FALSE)

# Upper and lower bounds
skinny_gibbs_genome_df_lymph <- 
  data.frame(dataset='lymph', n=n, dimension=dimension, no_chains=no_chains,
             W2L2UBmean_z=crn_bounds_z$W2L2_UBmean, W2L2UBsd_z=crn_bounds_z$W2L2_UBsd,
             W2L2LBmean_z=crn_bounds_z$W2L2_LBmean, W2L2LBsd_z=crn_bounds_z$W2L2_LBsd,
             W2L2UBmean_beta=crn_bounds_beta$W2L2_UBmean, W2L2UBsd_beta=crn_bounds_beta$W2L2_UBsd,
             W2L2LBmean_beta=crn_bounds_beta$W2L2_LBmean, W2L2LBsd_beta=crn_bounds_beta$W2L2_LBsd,
             W2L2UBs_z=I(list(crn_bounds_z$W2L2_UBs)), W2L2LBs_z=I(list(crn_bounds_z$W2L2_LBs)),
             W2L2UBs_beta=I(list(crn_bounds_beta$W2L2_UBs)), W2L2LBs_beta=I(list(crn_bounds_beta$W2L2_LBs)))
skinny_gibbs_genome_df <- rbind(skinny_gibbs_genome_df_malware, skinny_gibbs_genome_df_lymph)
# Saving data
# save(skinny_gibbs_genome_df, file="skinny_gibbs/skinny_gibbs_genome_df.RData")















# crn_exact_skinny_chain_sampler <- function(){
#   output <- exact_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
#   return(list('P'=output$exact_z, 'Q'=output$skinny_z))
# }
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
# 
# # Upper and lower bounds
# skinny_gibbs_genome_df <- 
#   data.frame(dataset='lymph', n=n, p=p, no_chains=no_chains,
#              coupled_W1hammingUBmean=coupled_W1hammingUBmean, coupled_W1hammingUBsd, coupled_W1hammingUBsd,
#              W1hammingLBmean=W1hammingLBmean, W1hammingLBsd=W1hammingLBsd,
#              indep_W1hammingUBmean=indep_W1hammingUBmean, indep_W1hammingUBsd=indep_W1hammingUBsd)
# 
# # Saving data
# # save(skinny_gibbs_genome_df, file="skinny_gibbs/skinny_gibbs_genome_df.RData")






# # if (!requireNamespace("BiocManager", quietly = TRUE))
# #   install.packages("BiocManager")
# # BiocManager::install("S4Vectors")
# # BiocManager::install("Biobase")
# # BiocManager::install("iBMQ")
# #load('/Users/niloybiswas/Downloads/E-GEOD-22575-atlasExperimentSummary.Rdata')
# #experimentSummary@listData$`A-AFFY-36`
# 
# 
# library("iBMQ")
# 
# data(PPA.liver)
# cutoff.liver <- calculateThreshold(PPA.liver, 0.2) 
# eqtl.liver <- eqtlFinder(PPA.liver, cutoff.liver) 
# hotspot.liver <- hotspotFinder(eqtl.liver,20)
# 
# data(phenotype.liver)
# phenotype.liver@assayData$exprs[5000,]
# dim(phenotype.liver@assayData$exprs)
# phenotype.liver@experimentData
# experimentData(phenotype.liver)
