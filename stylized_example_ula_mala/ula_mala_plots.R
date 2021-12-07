# Wasserstein's distance coupling plots
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

# Metric considered: L2
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}

# Analytic W2L2 bound of Durmus+Moulines
ula_w2_analytic_bound_bias <- function(cov_mat, std_mh){
  dimension <- dim(cov_mat)[1]
  eigenvalues <- eigen(cov_mat)$values
  L <- 1/min(eigenvalues)
  m <- 1/max(eigenvalues)
  k <- 2*m*L/(m+L)
  g <- std_mh^2/2
  w2_squared <- 2*k^(-1)*g^2*dimension*(2*L^2+L^4*(g^2/6 + g/m + k^(-1)*g + 4/(3*m)))
  # w2_squared <- 2*L^2*k^(-1)*g*(k^(-1)+g)*(2+L^2*g/m+L^2*(g^2)/6)*dimension
  return(sqrt(w2_squared))
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

# ULA target sampler
Sampler_ula_dist <- function(no_indep_samples, cov_mat, std_mh){
  dimension <- dim(cov_mat)[1]
  B <- diag(dimension)-chol2inv(chol(cov_mat))*(std_mh^2)/2
  ula_limit_cov_mat <- chol2inv(chol(diag(dimension)-cpp_prod(B,B)))*(std_mh^2)
  return(t(mvtnorm::rmvnorm(n = no_indep_samples, mean = rep(0,dimension), sigma = ula_limit_cov_mat)))
}

# p-Wasserstein considered: p=2
p = 2

################################################################################
# Plot: W2 bias estimate 
alpha <- 0.5
no_chains <- 10
burnin <- 1000
chain_length <- burnin+2e3
trajectory_bias_dimension_df_ula_mala <- data.frame()
for (dimension in seq(200,1000,200)){
  # Covariance matrix
  cov_mat <- cov_matrix(alpha, dimension)
  cov_mat_sqrt <- t(chol(cov_mat))
  inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)
  
  sigma_mh_ratio <- 0.5
  sigma_mh <- sigma_mh_ratio*(dimension^(-1/6))
  
  # Coupling-based bound
  coupled_chain_sampler_ula_mala <- 
    function(){
      coupled_chain <- 
        coupled_mala(matrix(rnorm(dimension), nrow = dimension, ncol = 1),
                     matrix(rnorm(dimension), nrow = dimension, ncol = 1), 
                     LogPdfP,LogPdfP, GradLogPdfP, GradLogPdfP, sigma_mh, 
                     chain_length, q_correction=FALSE, pi_correction=TRUE)
      return(list('P'=coupled_chain$P[-c(1:burnin),], 'Q'=coupled_chain$Q[-c(1:burnin),]))
    }
  crn_ula_mala <- wp_ub_estimate(coupled_chain_sampler_ula_mala, no_chains=no_chains,
                                 p=2, metric=metric_l2, parallel=TRUE)
  crn_coupling_ub_mean <- crn_ula_mala$wp_ub # mean(rowMeans(crn_ula_mala$wp_power_p_ub_tracjectories[,-c(1:burnin)]))
  crn_coupling_ub_sd <- crn_ula_mala$wp_ub_se # sd(rowMeans(crn_ula_mala$wp_power_p_ub_tracjectories[,-c(1:burnin)]))
  trajectory_bias_dimension_df_ula_mala <- 
    rbind(trajectory_bias_dimension_df_ula_mala, 
          data.frame(dimension=dimension, sigma_mh=sigma_mh, sigma_mh_ratio=sigma_mh_ratio,
                     metric_mean=crn_coupling_ub_mean, metric_sd=crn_coupling_ub_sd, 
                     type='averaged_ergodic'))
  
  # Analytic bound of Durmus and Moulines
  bias_ub_analytic_mean <- ula_w2_analytic_bound_bias(cov_mat, sigma_mh)
  bias_ub_analytic_sd <- 0
  trajectory_bias_dimension_df_ula_mala <- 
    rbind(trajectory_bias_dimension_df_ula_mala, 
          data.frame(dimension=dimension, sigma_mh=sigma_mh, sigma_mh_ratio=sigma_mh_ratio, 
                     metric_mean=bias_ub_analytic_mean, metric_sd=bias_ub_analytic_sd,
                     type='analytic'))
  
  # True W2L2 bias ULA
  true_w2_bias_mean <- ula_w2_bias(cov_mat, sigma_mh)
  true_w2_bias_sd <- 0
  trajectory_bias_dimension_df_ula_mala <- 
    rbind(trajectory_bias_dimension_df_ula_mala, 
          data.frame(dimension=dimension, sigma_mh=sigma_mh, sigma_mh_ratio=sigma_mh_ratio,
                     metric_mean=true_w2_bias_mean, metric_sd=true_w2_bias_sd,
                     type='true_w2'))
  
  # Independent coupling bias of ULA
  indep_w2_bias_mean <- indep_coupling_w2_bias(cov_mat, sigma_mh)
  indep_w2_bias_sd <- 0
  trajectory_bias_dimension_df_ula_mala <- 
    rbind(trajectory_bias_dimension_df_ula_mala, 
          data.frame(dimension=dimension, sigma_mh=sigma_mh, sigma_mh_ratio=sigma_mh_ratio,
                     metric_mean=indep_w2_bias_mean, metric_sd=indep_w2_bias_sd, 
                     type='indep_coupling'))
  
  print(c(dimension))
}

# save(trajectory_bias_dimension_df_ula_mala, file="stylized_example_ula_mala/trajectory_bias_dimension_comparison_df_ula_mala.RData")
# load("stylized_example_ula_mala/trajectory_bias_dimension_comparison_df_ula_mala.RData")

mvn_plot_dim_comparison <- 
  ggplot(trajectory_bias_dimension_df_ula_mala, 
         aes(x=dimension, y=metric_mean, linetype=type)) + 
  geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('$W_2$ upper bounds')) +
  geom_ribbon(aes(ymin=metric_mean-metric_sd, ymax=metric_mean+metric_sd),
              alpha=0.2, colour = NA) +
  scale_linetype_manual(name=TeX(''),
                        breaks = c("indep_coupling", "averaged_ergodic", "true_w2", "analytic"),
                        labels=unname(TeX(c('Indep. Coupling', 'CUB_2', 'True $W_2$', 'Durmus-Moulines'))),
                        values = c('dotted', 'solid', 'dashed', 'dotdash')) +
  # scale_x_continuous(breaks=seq(0,2e3,100), limits = c(0,150)) + 
  scale_y_continuous(limits = c(1e-2,1e2), breaks=10^seq(-2,3,1), trans='log10',
                     labels = scales::comma_format(accuracy = 0.01)) +
  # scale_y_log10() +
  theme_classic(base_size = 18) + 
  theme(legend.position = 'right', legend.key.width=unit(1.25,"cm")) +
  guides(linetype=guide_legend(ncol=1))
mvn_plot_dim_comparison
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/stylized_example_ula_mala/ula_mala_trajectory_dim_comparison_plot.pdf", plot = mvn_plot_dim_comparison, width = 8, height = 4)

