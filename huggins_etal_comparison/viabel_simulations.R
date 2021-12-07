rm(list = ls())
# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(transport)

library(doParallel)
registerDoParallel(cores = detectCores()-1)


setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source('kernels.R')
source('estimators.R')
source('rcpp_functions.R')
source('stylized_example_mvn/mvn_functions.R')

source('huggins_comparison/viabel_functions.R')


################################################################################
# W2L2 varying dimension
p <- 2 

alpha <- 0.5
no_chains <- 20
chain_length <- 1500
burnin <- 500
trajectory_dimension_df <- data.frame()
for (dimension in seq(2,10,2)){
  sigma_mh <- 0.5*(dimension^(-1/6))
  # Covariance matrix
  cov_mat <- cov_matrix(alpha, dimension)
  inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)
  
  coupled_chain_sampler_mala <- 
    function(){
      return(coupled_mala(SamplerQ(1),SamplerP(1), LogPdfQ,LogPdfP,
                          GradLogPdfQ, GradLogPdfP, sigma_mh, chain_length,
                          burnin))
    }
  crn_mala <- wp_ub_estimate(coupled_chain_sampler_mala, no_chains=no_chains,
                             p=2, metric=metric_l2, parallel=TRUE)
  
  crn_coupling_ub_mean <- crn_mala$wp_ub
  crn_coupling_ub_sd <- crn_mala$wp_ub_se
  
  # Linear program
  empirical_w2_power_p = rep(NA,no_chains)
  for (i in c(1:no_chains)){
    Q_samples <- t(SamplerQ(chain_length))
    Q_samples <- pp(Q_samples)
    P_samples <- t(SamplerP(chain_length))
    P_samples <- pp(P_samples)
    empirical_w2_power_p[i] <- transport::wasserstein(Q_samples, P_samples, p=p)^p
  }
  linear_program_mean = mean(empirical_w2_power_p)^(1/p)
  linear_program_sd = (sd(empirical_w2_power_p)/sqrt(no_chains))^(1/p)

  # Only for p=2
  indep_coupling_mean <- (2*dimension)^(1/p) 
  indep_coupling_sd <- 0
  
  # Only for p=2
  cov_mat_sqrt <- expm::sqrtm(cov_mat)
  true_W2_mean <- norm(cov_mat_sqrt-diag(dimension), type = 'F') # True W2L2 known
  true_W2_sd <- 0
  
  # Huggins bound
  no_samples <- chain_length*2
  unnormalised_target_logpdf <- LogPdfP
  LogPdfQ_normalized <- 
    function(x) {
      # dimension <- length(x)
      norm_constant <- -(dimension/2)*log(2*pi)
      return(norm_constant+LogPdfQ(x))
    }
  # NOTE: need proposal_logpdf to be normalized
  proposal_logpdf <- LogPdfQ_normalized 
  proposal_componentwise_variances <- diag(cov_mat)
  # elbo_proposal_logpdf <- LogPdfQ
  LogPdfP_normalized <-
    function(x) {
      # dimension <- length(x)
      norm_constant <- -(dimension/2)*log(2*pi)-(1/2)*log(det(cov_mat))
      return(norm_constant+LogPdfP(x))
    }
  elbo_proposal_logpdf <- LogPdfP_normalized
  
  huggins_W2L2_ubs <- rep(NA, no_chains)
  for(i in c(1:no_chains)){
    proposal_samples <- SamplerQ(no_samples)
    elbo_proposal_samples <- SamplerQ(no_samples)
    huggins_W2L2_ubs[i] <- 
      huggins_W2L2_ub_gaussian(unnormalised_target_logpdf, proposal_logpdf, proposal_samples, 
                               proposal_componentwise_variances, elbo_proposal_logpdf, elbo_proposal_samples)
  }
  huggins_mean = mean(huggins_W2L2_ubs^p)^(1/p)
  huggins_sd = (sd(huggins_W2L2_ubs^p)/sqrt(no_chains))^(1/p)
  
  trajectory_dimension_df <- 
    rbind(trajectory_dimension_df, 
          data.frame(metric_mean=true_W2_mean, metric_sd=true_W2_sd, type='true_W2L2', dimension=dimension),
          data.frame(metric_mean=crn_coupling_ub_mean, metric_sd=crn_coupling_ub_sd, type='mala_crn', dimension=dimension),
          data.frame(metric_mean=linear_program_mean, metric_sd=linear_program_sd, type='linear_program', dimension=dimension),
          data.frame(metric_mean=indep_coupling_mean, metric_sd=indep_coupling_sd, type='indep', dimension=dimension),
          data.frame(metric_mean=huggins_mean, metric_sd=huggins_sd, type='huggins_W2L2_ub', dimension=dimension))
  
  print(dimension)
}

mvn_plot_dim <- 
  ggplot(trajectory_dimension_df %>% filter(type != 'linear_program'), 
         aes(x=dimension, y=metric_mean, linetype=type)) + 
  geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('$W_2$ upper bounds')) +
  scale_linetype_manual(name=TeX(''), breaks = c('indep','mala_crn','huggins_W2L2_ub','true_W2L2'),
                        labels=unname(TeX(c('Indep. Coupling','CUB_2', 
                                            'Huggins et al.', 'True $W_2$'))),
                        values = c('dotted', 'solid', 'dotdash', 'dashed', 'longdash')) +
  geom_ribbon(aes(ymin=metric_mean-metric_sd/sqrt(no_chains),
                  ymax=metric_mean+metric_sd/sqrt(no_chains)), alpha=0.2, colour = NA) +
  # values = c(RColorBrewer::brewer.pal(5,name = 'Greys'))) +
  # scale_x_continuous(breaks=seq(0,2e3,200)) +
  # scale_y_continuous(limits = c(1,2100), trans = 'log10') +
  # scale_y_continuous(trans = 'log10') +
  # scale_y_continuous(breaks=seq(0,,1)) + 
  theme_classic(base_size = 18) + 
  theme(legend.key.width=unit(1.25,"cm"))
  # guides(linetype=guide_legend(nrow=2,byrow=TRUE))
mvn_plot_dim
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/comparison/huggins_comparison_plot.pdf", plot = mvn_plot_dim, width = 8, height = 4)






