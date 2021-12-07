# Functions to simulate bounds from the paper
# "Using Coupling Methods to Estimate Sample Quality of SDEs"
# https://epubs.siam.org/doi/abs/10.1137/20M1312009?journalCode=sjuqa3

# Distance between ULA and MALA targetting a gaussians

rm(list = ls())

# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(transport)
library(scales)

library(doParallel)
registerDoParallel(cores = detectCores()-1)

# Functions
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source('kernels.R')
source('estimators.R')
source('rcpp_functions.R')
source('stylized_example_mvn/mvn_functions.R')



# Metric considered: capped L2
metric_capped_l2 <- function(x,y){
  if(is.null(dim(x))){return(min(sum((x-y)^2)^0.5))} else {return(pmin(1,rowSums((x-y)^2)^0.5))}
}


# True W2L2 bias of ULA
ula_w2_bias <- function(cov_mat, std_mh){
  dimension <- dim(cov_mat)[1]
  B <- diag(dimension)-chol2inv(chol(cov_mat))*(std_mh^2)/2
  ula_limit_cov_mat <- chol2inv(chol(diag(dimension)-cpp_prod(B,B)))*(std_mh^2)
  sqrt_cov_mat <- expm::sqrtm(cov_mat)
  cross_term <- expm::sqrtm(sqrt_cov_mat%*%ula_limit_cov_mat%*%sqrt_cov_mat)
  return(sum(diag(cov_mat)+diag(ula_limit_cov_mat)-2*diag(cross_term))^0.5)
  # return(norm(ula_limit_cov_mat-cov_mat, type = 'F'))
}

# Indep coupling W2 bias of ULA
indep_coupling_w2_bias <- function(cov_mat, std_mh){
  dimension <- dim(cov_mat)[1]
  B <- diag(dimension)-chol2inv(chol(cov_mat))*(std_mh^2)/2
  ula_limit_cov_mat <- chol2inv(chol(diag(dimension)-cpp_prod(B,B)))*(std_mh^2)
  return(sqrt(sum(diag(ula_limit_cov_mat))+sum(diag(cov_mat))))
}

#### Dobson et al. paper functions ####
# Function to estimate contraction constant a
cntrctn_cnst_estimate <- function(omega_set,no_chains,burnin,chain_length){
  # Coupling-based bound
  coupled_chain_sampler_ula_ula <- 
    function(){
      init1 <- runif(dimension)*(omega_set[2]-omega_set[1])+omega_set[1]
      init2 <- runif(dimension)*(omega_set[2]-omega_set[1])+omega_set[1]
      # matrix(rnorm(dimension), nrow = dimension, ncol = 1)
      return(coupled_mala(matrix(init1, nrow = dimension, ncol = 1),
                          matrix(init2, nrow = dimension, ncol = 1), 
                          LogPdfP,LogPdfP,reflect_threshold=Inf,
                          GradLogPdfP, GradLogPdfP, sigma_mh, chain_length,
                          q_correction=FALSE, pi_correction=FALSE))
    }
  reflection_ula_ula <- wp_ub_estimate(coupled_chain_sampler_ula_ula, no_chains=no_chains,
                                 p=1, metric=metric_capped_l2, parallel=TRUE)
  
  meeting_time <- function(trajectory){
    if(sum(trajectory==0)>0){
      return(min(which(trajectory==0)))
    } else {
      return(length(trajectory))
    }
  }
  meeting_times <- 
    as.vector(apply(reflection_ula_ula$wp_power_p_ub_tracjectories, 1, meeting_time))
  
  minus_lambda_limit <- 
    sapply(c((1):chain_length), function(t) log(mean(meeting_times>t))/t)
  minus_lambda_limit <- minus_lambda_limit[minus_lambda_limit>-Inf]
  minus_lambda <- mean(minus_lambda_limit,trim=0.025)
  return(exp(minus_lambda*c((burnin+1):chain_length)))
}

# Function to estimate perturbation part
prtrb_bound_estimate <- function(no_chains, burnin, chain_length){
  crn_chain_sampler_ula_mala <- 
    function(){
      return(coupled_mala(SamplerP(1),
                          SamplerP(1), 
                          LogPdfP,LogPdfP,
                          GradLogPdfP, GradLogPdfP, sigma_mh, chain_length,
                          q_correction=FALSE, pi_correction=TRUE))
    }
  crn_ula_mala <- 
    wp_ub_estimate(crn_chain_sampler_ula_mala, no_chains=no_chains,
                   p=1, metric=metric_capped_l2, parallel=TRUE)
  
  means=apply(crn_ula_mala$wp_power_p_ub_tracjectories, 2, mean)[-c(1:burnin)]
  sds=apply(crn_ula_mala$wp_power_p_ub_tracjectories, 2, sd)[-c(1:burnin)]
  
  return(list('mean'=means, 'sd'=sds))
}

# Final bound
sample_quality_bound <- function(cntrctn_cnsts, prtrb_bounds, epsilon){
  
  return(list('mean'=(prtrb_bounds$mean+2*epsilon)/(1-cntrctn_cnsts), 
              'sd'=(prtrb_bounds$sd)/(1-cntrctn_cnsts)))
}


###############################################################################
# Plot 2: W1 bias estimate 
alpha <- 0.5
no_chains <- 100
burnin <- 1000
chain_length <- burnin+2e3

omega_set <- c(-4,4)
epsilon <- 1/no_chains

trajectory_bias_dimension_df_ula_mala <- data.frame()
for (dimension in seq(220,300,20)){
  # Covariance matrix
  cov_mat <- cov_matrix(alpha, dimension)
  cov_mat_sqrt <- t(chol(cov_mat))
  inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)
  
  sigma_mh_ratio <- 0.5
  sigma_mh <- sigma_mh_ratio*(dimension^(-1/6))
  
  # Coupling-based bound
  coupled_chain_sampler_ula_mala <- 
    function(){
      return(coupled_mala(matrix(rnorm(dimension), nrow = dimension, ncol = 1),
                          matrix(rnorm(dimension), nrow = dimension, ncol = 1), 
                          LogPdfP,LogPdfP,
                          GradLogPdfP, GradLogPdfP, sigma_mh, chain_length,
                          q_correction=FALSE, pi_correction=TRUE))
    }
  crn_ula_mala <- wp_ub_estimate(coupled_chain_sampler_ula_mala, no_chains=no_chains,
                                 p=1, metric=metric_capped_l2, parallel=TRUE)
  crn_coupling_ub_mean <- mean(rowMeans(crn_ula_mala$wp_power_p_ub_tracjectories[,-c(1:burnin)]))
  crn_coupling_ub_sd <- sd(rowMeans(crn_ula_mala$wp_power_p_ub_tracjectories[,-c(1:burnin)]))
  trajectory_bias_dimension_df_ula_mala <- 
    rbind(trajectory_bias_dimension_df_ula_mala, 
          data.frame(dimension=dimension, sigma_mh=sigma_mh, sigma_mh_ratio=sigma_mh_ratio,
                     metric_mean=crn_coupling_ub_mean, metric_sd=crn_coupling_ub_sd, 
                     type='averaged_ergodic'))
  
  # Dobson bound
  cntrctn_cnsts <- cntrctn_cnst_estimate(omega_set,no_chains,burnin,chain_length)
  prtrb_bounds <- prtrb_bound_estimate(no_chains,burnin,chain_length)
  dobson_ub <- sample_quality_bound(cntrctn_cnsts, prtrb_bounds, epsilon)
  
  trajectory_bias_dimension_df_ula_mala <-
    rbind(trajectory_bias_dimension_df_ula_mala,
          data.frame(dimension=dimension, sigma_mh=sigma_mh, sigma_mh_ratio=sigma_mh_ratio,
                     metric_mean=mean(dobson_ub$mean), metric_sd=sd(dobson_ub$sd),
                     type='dobson'))
  
  print(c(dimension))
}

# save(trajectory_bias_dimension_df_ula_mala, file="dobson_paper_comparison/dobson_comparison_df_ula_mala.RData")
# load("dobson_paper_comparison/dobson_comparison_df_ula_mala.RData")


# Plot
mvn_plot_dim_comparison <- 
  ggplot(trajectory_bias_dimension_df_ula_mala, 
         aes(x=dimension, y=metric_mean, linetype=type)) + 
  geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('$W_1$ upper bounds')) +
  # geom_ribbon(aes(ymin=metric_mean-metric_sd, ymax=metric_mean+metric_sd), alpha=0.2, colour = NA) +
  scale_linetype_manual(name=TeX(''),
                        breaks = c("averaged_ergodic", "dobson"),
                        labels=unname(TeX(c('CUB_1', 'Dobson et al.'))),
                        values = c('solid', 'dotdash')) +
  # scale_x_continuous(breaks=seq(0,2e3,100), limits = c(0,150)) + 
  # scale_y_continuous(limits = c(1e-3,25e2), breaks=10^seq(-2,3,1), trans='log10',
  #                    labels = scales::comma_format(accuracy = 0.01)) +
  # scale_y_log10() +
  # scale_y_continuous(limits = c(0.5,1.25)) + 
  theme_classic(base_size = 18) + 
  theme(legend.position = 'right', legend.key.width=unit(1.25,"cm")) +
  guides(linetype=guide_legend(ncol=1))
mvn_plot_dim_comparison
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/comparison/dobson_comparison_plot.pdf", plot = mvn_plot_dim_comparison, width = 8, height = 4)




