# Wasserstein's distance coupling plots
# Distance between two Gaussians

rm(list = ls())

# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(transport)

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

# p-Wasserstein considered: p=2 (don not change)
p = 2

################################################################################
# Plot 1: dimension=100
dimension <- 100
alpha <- 0.5

# Covariance matrix
cov_mat <- cov_matrix(alpha, dimension)
inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)

# Coupled MALA
chain_length <- 1000
sigma_mh <- 0.5*(dimension^(-1/6))
coupled_chain_sampler_mala <- 
  function(){
    return(coupled_mala(SamplerQ(1),SamplerP(1), LogPdfQ,LogPdfP,
                 GradLogPdfQ, GradLogPdfP, sigma_mh, chain_length))
  }
no_chains <- 10
crn_mala <- wp_ub_estimate(coupled_chain_sampler_mala, no_chains=no_chains,
                           p=p, metric=metric_l2, parallel=FALSE)
rownames(crn_mala$wp_power_p_ub_tracjectories) <- c(1:no_chains)

power_p_ergodic_means <- apply(crn_mala$wp_power_p_ub_tracjectories,1,cummean)

single_ergodic_trajectory_df <- 
  data.frame(t=c(1:chain_length), power_p_ergodic_means) %>%
  reshape2::melt(id.vars = c("t"), value.name = c('metric_power_p')) %>%
  dplyr::rename(trajectory=variable) %>%
  dplyr::mutate(type='single_ergodic', trajectory=as.integer(trajectory))

# Trajectory means and sds
d_coupling <- single_ergodic_trajectory_df %>% 
  dplyr::filter(type=='single_ergodic') %>%
  dplyr::group_by(t) %>%
  dplyr::summarise(metric_mean=mean(metric_power_p)^(1/p), metric_sd=(sd(metric_power_p)/sqrt(no_chains))^(1/p))
trajectory_df <- 
  data.frame(t=c(1:chain_length), metric_mean=d_coupling$metric_mean,
             metric_sd=d_coupling$metric_sd, type='averaged_ergodic', 
             no_chains=no_chains)



# Linear program
empirical_w2_power_p = rep(NA,no_chains)
for (i in c(1:no_chains)){
  Q_samples <- t(SamplerQ(chain_length))
  Q_samples <- pp(Q_samples)
  P_samples <- t(SamplerP(chain_length))
  P_samples <- pp(P_samples)
  empirical_w2_power_p[i] <- transport::wasserstein(Q_samples, P_samples, p=p)^p
}
trajectory_df <- 
  rbind(trajectory_df, 
        data.frame(t=c(1:chain_length), metric_mean=(mean(empirical_w2_power_p))^(1/p),
                   metric_sd=(sd(empirical_w2_power_p)/sqrt(no_chains))^(1/p), type='linear_program',
                   no_chains=no_chains))

# Independent coupling
# indep_coupling_samples <- metric_l2(t(SamplerQ(chain_length*no_chains)), t(SamplerP(chain_length*no_chains)))
# indep_coupling_mean <- mean(indep_coupling_samples)
# indep_coupling_sd <- sd(indep_coupling_samples)/sqrt(chain_length*no_chains)
indep_coupling_mean <- (2*dimension)^(1/2)
indep_coupling_sd <- 0

# indep_coupling_mean <- # As we initialize the chain from the indep coupling
#   (trajectory_df %>% filter(t==1, type=='averaged_ergodic'))$metric_mean
# indep_coupling_sd <- # As we initialize the chain from the indep coupling
#   (trajectory_df %>% filter(t==1, type=='averaged_ergodic'))$metric_sd
trajectory_df <- 
  rbind(trajectory_df, 
        data.frame(t=c(1:chain_length), metric_mean=indep_coupling_mean,
                   metric_sd=indep_coupling_sd, type='indep_coupling',
                   no_chains=no_chains))

# True W2L2
cov_mat_sqrt <- expm::sqrtm(cov_mat)
true_W2_squared_mean <- norm(cov_mat_sqrt-diag(dimension), type = 'F') # True W2L2 known
true_W2_squared_sd <- 0
trajectory_df <- 
  rbind(trajectory_df, 
        data.frame(t=c(1:chain_length), metric_mean=true_W2_squared_mean,
                   metric_sd=true_W2_squared_sd, type='true_w2',
                   no_chains=no_chains))

# save(trajectory_df, file="stylized_example_mvn/trajectory_df.RData")
# load("stylized_example_mvn/trajectory_df.RData")

mvn_plot <- 
  ggplot(trajectory_df, aes(x=t, y=metric_mean, linetype=type)) + 
  geom_line(size=1) + xlab(TeX('Trajectory length $T$')) + ylab(TeX('$W_2$ upper bounds')) +
  geom_ribbon(aes(ymin=metric_mean-metric_sd, ymax=metric_mean+metric_sd), 
              alpha=0.2, colour = NA) +
  scale_linetype_manual(name=TeX(''),
                        breaks = c("indep_coupling", "averaged_ergodic", "linear_program", "true_w2"),
                        labels=unname(TeX(c('Indep. Coupling', 'CUB_2', 'Empirical W_2', 'True $W_2$'))),
                        values = c('dotted', 'solid', 'dotdash', 'dashed')) +
  scale_y_continuous(limits = c(0,18)) +
  theme_classic(base_size = 12) + 
  theme(legend.position = 'bottom', legend.key.width=unit(1.25,"cm")) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))
mvn_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/stylized_example_mvn/mvn_trajectory_plot.pdf", plot = mvn_plot, width = 4, height = 4)


################################################################################
# Plot 2: W2 varying dimension
alpha <- 0.5
no_chains <- 5
chain_length <- 1000
burnin <- 0
trajectory_dimension_df <- data.frame()
for (dimension in seq(200,1000,200)){
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
                             p=2, metric=metric_l2, parallel=FALSE)
  
  crn_coupling_ub_mean <- crn_mala$wp_ub # mean(rowMeans(crn_mala$wp_power_p_ub_tracjectories))
  crn_coupling_ub_sd <- crn_mala$wp_ub_se # sd(rowMeans(crn_mala$wp_power_p_ub_tracjectories))
  
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
  
  indep_coupling_mean <- (2*dimension)^(1/p)
  indep_coupling_sd <- 0
  
  cov_mat_sqrt <- expm::sqrtm(cov_mat)
  true_W2_mean <- norm(cov_mat_sqrt-diag(dimension), type = 'F') # True W2L2 known
  true_W2_sd <- 0
  
  trajectory_dimension_df <- 
    rbind(trajectory_dimension_df, 
          data.frame(metric_mean=true_W2_mean, metric_sd=true_W2_sd, type='true_W2L2', dimension=dimension),
          data.frame(metric_mean=crn_coupling_ub_mean, metric_sd=crn_coupling_ub_sd, type='mala_crn', dimension=dimension),
          data.frame(metric_mean=linear_program_mean, metric_sd=linear_program_sd, type='linear_program', dimension=dimension),
          data.frame(metric_mean=indep_coupling_mean, metric_sd=indep_coupling_sd, type='indep', dimension=dimension))
  print(dimension)
}

# save(trajectory_dimension_df, file="stylized_example_mvn/trajectory_dimension_df.RData")
# load("stylized_example_mvn/trajectory_dimension_df.RData")

mvn_plot_dim <- 
  ggplot(trajectory_dimension_df, aes(x=dimension, y=metric_mean, linetype=type)) + 
  geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('$W_2$ upper bounds')) +
  scale_linetype_manual(name=TeX(''), breaks = c('indep','mala_crn', 'linear_program','true_W2L2'),
                        labels=unname(TeX(c('Indep. Coupling','CUB_2', 
                                            'Empirical $W_2$', 'True $W_2$'))),
                        values = c('dotted', 'solid', 'dotdash', 'dashed')) +
  geom_ribbon(aes(ymin=metric_mean-metric_sd,
                  ymax=metric_mean+metric_sd), alpha=0.2, colour = NA) +
  # values = c(RColorBrewer::brewer.pal(5,name = 'Greys'))) +
  scale_x_continuous(breaks=seq(0,2e3,200)) +
  scale_y_continuous(limits = c(0,50)) +
  # scale_y_continuous(breaks=seq(0,,1)) + 
  theme_classic(base_size = 12) + 
  theme(legend.position = 'bottom', legend.key.width=unit(1.25,"cm")) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))
mvn_plot_dim
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/stylized_example_mvn/mvn_trajectory_dim_plot.pdf", plot = mvn_plot_dim, width = 4, height = 4)








