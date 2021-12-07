# Sinkhorn distance comparison
rm(list = ls())

# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(transport)
library(Barycenter)


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

p <- 2

################################################################################
# Plot 1: dimension=10
dimension <- 10
alpha <- 0.5

# Covariance matrix
cov_mat <- cov_matrix(alpha, dimension)
inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)
cov_mat_sqrt <- t(chol(cov_mat))

# Coupled MALA
chain_length <- 500
sigma_mh <- 0.5
coupled_chain_sampler_mala <- 
  function(){
    return(coupled_mala(SamplerQ(1),SamplerP(1), LogPdfQ,LogPdfP,
                        GradLogPdfQ, GradLogPdfP, sigma_mh, chain_length))
  }
no_chains <- 10

start.time <- proc.time()[[3]]
crn_mala <- wp_ub_estimate(coupled_chain_sampler_mala, no_chains=no_chains,
                           p=p, metric=metric_l2, parallel=TRUE)
end.time <- proc.time()[[3]]
time_taken_coupling <- as.vector(round(end.time - start.time,2))

### Sinkhorn with independent samples
n_samples <- chain_length*no_chains
q_samples <- SamplerQ(n_samples)
pi_samples <- SamplerP(n_samples)

cost_mat <- matrix(NA, n_samples, n_samples)
for (i in 1:(n_samples)){
  for (j in 1:(n_samples)){
    cost_mat[i,j] <- metric_l2(q_samples[,i],pi_samples[,j])^2
  }
  print(i)
}

nrepeats <- 1
lambda_grid <- seq(0.25,10,0.25)
sinkhorn_distances <- rep(NA, length(lambda_grid))
time_taken_sinkhorn <- rep(NA, length(lambda_grid))
for (i in 1:length(lambda_grid)){
  time_taken_repeats <- rep(NA, nrepeats)
  sinkhorn_repeats <- rep(NA, nrepeats)
  for(j in 1:nrepeats){
    start.time <- proc.time()[[3]]
    output <- Sinkhorn(matrix(1/n_samples,n_samples,1), matrix(1/n_samples,n_samples,1), 
                       cost_mat, lambda=lambda_grid[i])
    end.time <- proc.time()[[3]]
    time_taken_repeats[j] <- round(end.time - start.time,2)
    sinkhorn_repeats[j] <- (output$Distance)^0.5
  }
  time_taken_sinkhorn[i] <- mean(time_taken_repeats)
  sinkhorn_distances[i] <- mean(sinkhorn_repeats)
  print(i)
}

### Exact Wasserstein solver
Q_samples <- t(q_samples)
Q_samples <- pp(Q_samples)
P_samples <- t(pi_samples)
P_samples <- pp(P_samples)
start.time <- proc.time()[[3]]
wass_solver <- transport::wasserstein(Q_samples, P_samples, p=p)
end.time <- proc.time()[[3]]
time_taken_wass_solver <- as.vector(round(end.time - start.time,2))

indep_coupling_distance <- mean(metric_l2(t(q_samples), t(pi_samples)))
true_W2 <- norm(cov_mat_sqrt-diag(dimension), type = 'F') # True W2L2 known

burnin <- 100
coupling_based_W2ub <- mean(crn_mala$wp_power_p_ub_mean_tracjectory[burnin:chain_length])^(1/p)

sinkhorn_df <- 
  rbind(data.frame(lambda=lambda_grid, metric=sinkhorn_distances, time=time_taken_sinkhorn, algo='sinkhorn'),
        data.frame(lambda=lambda_grid, metric=wass_solver, time=time_taken_wass_solver, algo='wass_solver'),
        data.frame(lambda=lambda_grid, metric=coupling_based_W2ub, time=time_taken_coupling, algo='mala_crn'),
        data.frame(lambda=lambda_grid, metric=indep_coupling_distance, time=NA, algo='indep'),
        data.frame(lambda=lambda_grid, metric=true_W2, time=NA, algo='trueW2L2'))
# save(sinkhorn_df, file="sinkhorn_comparison/sinkhorn_df.RData")
# load('sinkhorn_comparison/sinkhorn_df.RData')

# Distance Plot
plot.breaks <- c('sinkhorn','wass_solver', 'mala_crn','indep','trueW2L2')
plot.labels <- unname(TeX(c('Sinkhorn','Wasserstein distance solver','MALA CRN','Indep. Coupling','$W_2(P,Q)$')))
sinkhorn_plot <- ggplot(sinkhorn_df %>% filter(algo != 'wass_solver'), 
                        aes(x=lambda, y=metric, linetype=algo, color=algo)) + 
  geom_line(size=1) + xlab(TeX('Entropic regularisation parameter $\\lambda$')) + ylab(TeX('Distance')) +
  scale_linetype_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels,
                        values = c('solid','dotted','solid', 'dotted', 'dotdash')) +
  scale_color_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels,
                     values = c('gray','gray','black', 'black', 'black')) +
  # values = c(RColorBrewer::brewer.pal(5,name = 'Greys'))) +
  # scale_x_continuous(breaks=seq(0,2e3,200)) +
  scale_y_continuous(breaks=seq(0,5,1), limits = c(0,5)) + 
  theme_classic(base_size = 12) + 
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(linetype=guide_legend(nrow=2,byrow=TRUE))
sinkhorn_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/sinkhorn/sinkhorn_plot.pdf", plot = sinkhorn_plot, width = 4, height = 4)


# Time Plot
plot.breaks <- c('sinkhorn','mala_crn','indep','trueW2L2')
plot.labels <- unname(TeX(c('Sinkhorn','MALA CRN','Indep. Coupling','$W_2(P,Q)$')))
sinkhorn_time_plot <- 
  ggplot(sinkhorn_df %>% filter(algo %in% c("sinkhorn","mala_crn")), 
         aes(x=lambda, y=time, color=algo)) + geom_line(size=1) +
  xlab(TeX('Entropic regularisation parameter $\\lambda$')) + 
  ylab(TeX('Time (seconds)')) +
  # scale_y_continuous(breaks=seq(0,1000,200), limits = c(0,600)) +
  scale_y_log10() +
  scale_color_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels, values = c('gray','black')) +
  theme_classic(base_size = 12) + 
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
sinkhorn_time_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/sinkhorn/sinkhorn_times_plot.pdf", plot = sinkhorn_time_plot, width = 4, height = 4)


################################################################################
# Extra plots
################################################################################
# # Plot 2: dimension=100
# dimension <- 100
# alpha <- 0.5
# 
# # Covariance matrix
# cov_mat <- cov_matrix(alpha, dimension)
# inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)
# cov_mat_sqrt <- t(chol(cov_mat))
# 
# # Coupled MALA
# chain_length <- 2000
# sigma_mh <- 0.5
# coupled_chain_sampler_mala <- 
#   function(){
#     return(coupled_mala(SamplerP(1),SamplerP(1), LogPdfP,LogPdfP,
#                         GradLogPdfP, GradLogPdfP, sigma_mh, chain_length))
#   }
# no_chains <- 1
# 
# start.time <- proc.time()[[3]]
# crn_mala <- wp_ub_estimate(coupled_chain_sampler_mala, no_chains=no_chains,
#                            p=p, metric=metric_l2, parallel=TRUE)
# end.time <- proc.time()[[3]]
# time_taken_coupling <- as.vector(round(end.time - start.time,2))
# 
# ### Sinkhorn with independent samples
# n_samples <- chain_length*no_chains-burnin
# q_samples <- SamplerP(n_samples)
# pi_samples <- SamplerP(n_samples)
# 
# cost_mat <- matrix(NA, n_samples, n_samples)
# for (i in 1:(n_samples)){
#   for (j in 1:(n_samples)){
#     cost_mat[i,j] <- metric_l2(q_samples[,i],pi_samples[,j])^2
#   }
#   print(i)
# }
# 
# nrepeats <- 1
# lambda_grid <- c(seq(0.25, 1, 0.25),seq(2,30,1))
# sinkhorn_distances <- rep(NA, length(lambda_grid))
# time_taken_sinkhorn <- rep(NA, length(lambda_grid))
# for (i in 1:length(lambda_grid)){
#   time_taken_repeats <- rep(NA, nrepeats)
#   sinkhorn_repeats <- rep(NA, nrepeats)
#   for(j in 1:nrepeats){
#     start.time <- proc.time()[[3]]
#     output <- Sinkhorn(matrix(1/n_samples,n_samples,1), matrix(1/n_samples,n_samples,1), 
#                        cost_mat, lambda=lambda_grid[i])
#     end.time <- proc.time()[[3]]
#     time_taken_repeats[j] <- round(end.time - start.time,2)
#     sinkhorn_repeats[j] <- (output$Distance)^0.5
#   }
#   time_taken_sinkhorn[i] <- mean(time_taken_repeats)
#   sinkhorn_distances[i] <- mean(sinkhorn_repeats)
#   print(i)
# }
# 
# indep_coupling_distance <- mean(metric_l2(t(q_samples), t(pi_samples)))
# true_W2 <- 0 # True W2L2 zero by construction
# burnin <- 1e3
# coupling_based_W2ub <- mean(crn_mala$wp_power_p_ub_mean_tracjectory[burnin:chain_length])^(1/p)
# sinkhorn_df_zero <- 
#   rbind(data.frame(lambda=lambda_grid, metric=sinkhorn_distances, time=time_taken_sinkhorn, algo='sinkhorn'),
#         data.frame(lambda=lambda_grid, metric=coupling_based_W2ub, time=time_taken_coupling, algo='mala_crn'),
#         data.frame(lambda=lambda_grid, metric=indep_coupling_distance, time=NA, algo='indep'),
#         data.frame(lambda=lambda_grid, metric=true_W2, time=NA, algo='trueW2L2'))
# # save(sinkhorn_df_zero, file="sinkhorn_comparison/sinkhorn_df_zero.RData")
# # load('sinkhorn_comparison/sinkhorn_df_zero.RData')
# 
# # Distance Plot
# plot.breaks <- c('sinkhorn','mala_crn','indep', "trueW2L2")
# plot.labels <- unname(TeX(c('Sinkhorn','MALA CRN','Indep. Coupling','$W_2(P,Q)$ Optimal')))
# sinkhorn_zero_plot <- 
#   ggplot(sinkhorn_df_zero %>% filter(algo != "trueW2L2"), 
#          aes(x=lambda, y=metric, linetype=algo, color=algo)) + 
#   geom_line(size=1) + xlab(TeX('Entropic regularisation parameter $\\lambda$')) + ylab(TeX('Distance')) +
#   scale_linetype_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels,
#                         values = c('solid','solid', 'dotted', 'dotdash')) +
#   scale_color_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels,
#                      values = c('gray','black', 'black', 'black')) +
#   # values = c(RColorBrewer::brewer.pal(5,name = 'Greys'))) +
#   scale_x_continuous(breaks=seq(0,2e3,5), limits = c(0,25)) +
#   scale_y_continuous(breaks=seq(0,20,5), limits = c(0,15)) + 
#   theme_classic(base_size = 12) + 
#   theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
#   guides(linetype=guide_legend(nrow=2,byrow=TRUE))
# sinkhorn_zero_plot
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/sinkhorn/sinkhorn_zero_plot.pdf", plot = sinkhorn_zero_plot, width = 4, height = 4)
# 
# 
# 
# ################################################################################
# # Plot 2: varying dimension
# alpha <- 0.5
# 
# sinkhorn_df_varying_dim <- data.frame()
# for (dimension in seq(10,50,10)){
#   # Covariance matrix
#   cov_mat <- cov_matrix(alpha, dimension)
#   inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)
#   cov_mat_sqrt <- t(chol(cov_mat))
#   
#   chain_length <- 2000
#   sigma_mh <- 0.5*(dimension)^(-1/6)
#   coupled_chain_sampler_mala <- 
#     function(){
#       return(coupled_mala(SamplerP(1),SamplerP(1), LogPdfP,LogPdfP,
#                           GradLogPdfP, GradLogPdfP, sigma_mh, chain_length))
#     }
#   no_chains <- 1
#   crn_mala <- wp_ub_estimate(coupled_chain_sampler_mala, no_chains=no_chains,
#                              p=p, metric=metric_l2, parallel=TRUE)
#   
#   start.time <- proc.time()[[3]]
#   crn_mala <- wp_ub_estimate(coupled_chain_sampler_mala, no_chains=no_chains,
#                              p=p, metric=metric_l2, parallel=TRUE)
#   end.time <- proc.time()[[3]]
#   time_taken_coupling <- as.vector(round(end.time - start.time,2))
#   
#   burnin <- 1e3
#   
#   # Independent samples
#   n_samples <- chain_length*no_chains #-burnin
#   q_samples <- SamplerP(n_samples)
#   pi_samples <- SamplerP(n_samples)
#   
#   cost_mat <- matrix(NA, n_samples, n_samples)
#   for (i in 1:(n_samples)){
#     for (j in 1:(n_samples)){
#       cost_mat[i,j] <- metric_l2(q_samples[,i],pi_samples[,j])^2
#     }
#     print(i)
#   }
#   
#   ### Sinkhorn with independent samples
#   nrepeats <- 1
#   lambda_grid <- c(0.25)
#   sinkhorn_distances <- rep(NA, length(lambda_grid))
#   time_taken_sinkhorn <- rep(NA, length(lambda_grid))
#   for (i in 1:length(lambda_grid)){
#     time_taken_repeats <- rep(NA, nrepeats)
#     sinkhorn_repeats <- rep(NA, nrepeats)
#     for(j in 1:nrepeats){
#       start.time <- proc.time()[[3]]
#       output <- Sinkhorn(matrix(1/n_samples,n_samples,1), matrix(1/n_samples,n_samples,1), 
#                          cost_mat, lambda=lambda_grid[i])
#       end.time <- proc.time()[[3]]
#       time_taken_repeats[j] <- round(end.time - start.time,2)
#       sinkhorn_repeats[j] <- (output$Distance)^0.5
#     }
#     time_taken_sinkhorn[i] <- mean(time_taken_repeats)
#     sinkhorn_distances[i] <- mean(sinkhorn_repeats)
#     print(i)
#   }
#   
#   ### Exact Wasserstein solver
#   Q_samples <- t(q_samples)
#   Q_samples <- pp(Q_samples)
#   P_samples <- t(pi_samples)
#   P_samples <- pp(P_samples)
#   start.time <- proc.time()[[3]]
#   wass_solver <- transport::wasserstein(Q_samples, P_samples, p=p)
#   end.time <- proc.time()[[3]]
#   time_taken_wass_solver <- as.vector(round(end.time - start.time,2))
#   
#   indep_coupling_distance <- mean(metric_l2(t(q_samples), t(pi_samples)))
#   true_W2 <- 0 # True W2L2 zero by construction
#   coupling_based_W2ub <- mean(crn_mala$wp_power_p_ub_mean_tracjectory[burnin:chain_length])^(1/p)
#   sinkhorn_df_varying_dim <- 
#     rbind(sinkhorn_df_varying_dim,
#           data.frame(lambda=lambda_grid, metric=sinkhorn_distances, time=time_taken_sinkhorn, algo='sinkhorn', dimension=dimension),
#           data.frame(lambda=lambda_grid, metric=wass_solver, time=time_taken_wass_solver, algo='wass_solver', dimension=dimension),
#           data.frame(lambda=lambda_grid, metric=coupling_based_W2ub, time=time_taken_coupling, algo='mala_crn', dimension=dimension),
#           data.frame(lambda=lambda_grid, metric=indep_coupling_distance, time=NA, algo='indep', dimension=dimension),
#           data.frame(lambda=lambda_grid, metric=true_W2, time=NA, algo='trueW2L2', dimension=dimension))
#   
#   print(dimension)
# }
# 
# # save(sinkhorn_df_varying_dim, file="sinkhorn_comparison/sinkhorn_df_varying_dim.RData")
# # load('sinkhorn_comparison/sinkhorn_df_varying_dim.RData')
# 
# 
# # Distance Plot
# plot.breaks <- c("wass_solver",'mala_crn','indep')
# plot.labels <- 
#   unname(TeX(c('$W_2(P_N,Q_N)$','MALA CRN','Indep. Coupling')))
# sinkhorn_dim_plot <- 
#   ggplot(sinkhorn_df_varying_dim %>% filter(algo != "trueW2L2" & algo != "sinkhorn"), 
#          aes(x=dimension, y=metric, linetype=algo, color=algo)) + 
#   geom_line(size=1) + xlab(TeX('Dimension $d$')) + ylab(TeX('Distance')) +
#   scale_linetype_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels,
#                         values = c('solid','solid', 'dotted')) +
#   scale_color_manual(name=TeX('Algorithm'), breaks = plot.breaks, labels=plot.labels,
#                      values = c('gray','black', 'black')) +
#   # values = c(RColorBrewer::brewer.pal(5,name = 'Greys'))) +
#   # scale_x_continuous(breaks=seq(0,2e3,5), limits = c(0,25)) +
#   scale_y_continuous(breaks=seq(0,20,5), limits = c(0,10)) + 
#   theme_classic(base_size = 12) + 
#   theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
#   guides(linetype=guide_legend(nrow=2,byrow=TRUE))
# sinkhorn_dim_plot
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/sinkhorn/sinkhorn_dim_plot.pdf", plot = sinkhorn_dim_plot, width = 4, height = 4)


