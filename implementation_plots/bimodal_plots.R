rm(list = ls())
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source(file = 'kernels.R')
source(file = 'estimators.R')
source(file = 'implementation_plots/mvn_mixture_functions.R')

# Libraries
library(dplyr)
library(ggplot2)
library(latex2exp)
library(mcmcse)

# library(dplyr)
library(doParallel)
registerDoParallel(cores = detectCores()-1)

# Metric considered: L2
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}

# Effective sample size calculator
ess_calculator <- function(data){
  if(is.null(dim(data))){
    if(length(data)<100) {return(NA)}
    return(mcmcse::ess(data))
  } else {
    if(dim(data)[1]<100) {return(NA)}
    return(mcmcse::ess(data))
  }
  return(mcmcse::ess(data))
}

################################################################################
#  Targets: mixture of gaussians
dimension <- 4
no_components_P <- 2
weights_P <- rep(1,no_components_P)
weights_P <- weights_P/sum(weights_P)
means_P <- cbind(rep(-1,dimension),rep(1,dimension))

weights_Q <- rep(1,no_components_P)
weights_Q <- weights_Q/sum(weights_Q)
means_Q <- cbind(rep(1,dimension),rep(1,dimension))

# # Marginal MALA chains
# chain_length <- 1000
# init_state_P <- SamplerP(1)[1,]
# sigma_mh_P <- 1.5
# mala_output_P <- mala(init_state_P, LogPdfP, GradLogPdfP, sigma_mh_P, chain_length)
# mean(mala_output_P[2:chain_length,1]!=mala_output_P[1:(chain_length-1),1])
# init_state_Q <- SamplerQ(1)[1,]
# sigma_mh_Q <- 1.5
# mala_output_Q <- mala(init_state_Q, LogPdfQ, GradLogPdfQ, sigma_mh_Q, chain_length)
# mean(mala_output_Q[2:chain_length,1]!=mala_output_Q[1:(chain_length-1),1])
# 
# par(mfrow=c(2,2))
# matplot(mala_output_P[,1], type = 'l')
# hist(mala_output_P[,1])
# matplot(mala_output_Q[,1], type = 'l')
# hist(mala_output_Q[,1])

# CRN Coupled MALA
chain_length <- 3000
sigma_mh_P <- 1*(dimension^(-1/6))
sigma_mh_Q <- 1*(dimension^(-1/6))
sigma_mh <- c(sigma_mh_P,sigma_mh_Q)

coupled_crn_mala_sampler <- 
  function(){
    # return(coupled_mala(SamplerP(1)[1,], SamplerQ(1)[1,], LogPdfP, LogPdfQ,
    #                     GradLogPdfP, GradLogPdfQ, sigma_mh, chain_length))
    return(coupled_mala(rep(1,dimension), rep(1,dimension), LogPdfP, LogPdfQ,
                        GradLogPdfP, GradLogPdfQ, sigma_mh, chain_length))
  }
p <- 1
crn_ula_mala_1_chain <- 
  wp_ub_estimate(coupled_crn_mala_sampler, no_chains=1,
                 p=p, metric=metric_l2, parallel=TRUE)
single_chain_ergodic_ess <- 
  sapply(c(1:chain_length), function(x) 
    ess_calculator(((crn_ula_mala_1_chain$wp_power_p_ub_mean_tracjectory)^(1/p))[1:x]))
single_chain_ergodic_sd <- 
  sapply(c(1:chain_length), function(x) 
    sd(((crn_ula_mala_1_chain$wp_power_p_ub_mean_tracjectory)^(1/p))[1:x]))

no_chains <- 100
crn_ula_mala_multi_chain_trajectory <- 
  wp_ub_estimate(coupled_crn_mala_sampler, no_chains=no_chains,
                 p=p, metric=metric_l2, parallel=TRUE)
multi_chain_sd <- 
  apply((crn_ula_mala_multi_chain_trajectory$wp_power_p_ub_tracjectories)^(1/p), 2, sd)

bimodal_trajectory_df <-
  rbind(data.frame(t=c(1:chain_length), type='1_chain', trajectory=1,
             metric=(crn_ula_mala_1_chain$wp_power_p_ub_mean_tracjectory)^(1/p),
             metric_sd=single_chain_ergodic_ess/single_chain_ergodic_sd),
  data.frame(t=c(1:chain_length), type='100_chains', trajectory=2,
             metric=(crn_ula_mala_multi_chain_trajectory$wp_power_p_ub_mean_tracjectory)^(1/p),
             metric_sd=apply((crn_ula_mala_multi_chain_trajectory$wp_power_p_ub_tracjectories)^(1/p), 2, sd)/sqrt(no_chains)),
               #sd((crn_ula_mala_multi_chain_trajectory$wp_power_p_ub_mean_tracjectory)^(1/p))),
  data.frame(t=c(1:chain_length), type='1_chain_ergodic', trajectory=3,
             metric=cummean((crn_ula_mala_1_chain$wp_power_p_ub_mean_tracjectory)^(1/p)),
             metric_sd=multi_chain_sd/sqrt(no_chains)))

# True W1L1 distance
Q_samples <- SamplerQ(15000)
Q_samples <- transport::pp(Q_samples)
P_samples <- SamplerP(15000)
P_samples <- transport::pp(P_samples)
# true_w1 <- transport::wasserstein(Q_samples, P_samples)
true_w1 <- 2.027868 # known value obtained from transport::wasserstein(Q_samples, P_samples)
bimodal_trajectory_df <- 
  rbind(bimodal_trajectory_df, 
        data.frame(t=c(1:chain_length), type='true_w1', trajectory=4, 
                   metric=true_w1, metric_sd=NA))

# Independent coupling
indep_coupling <- mean(metric_l2(SamplerQ(1000), SamplerP(1000)))
# indep_coupling <- # As we initialize the chain from the indep coupling
#   (bimodal_trajectory_df %>% filter(t==1, type == '100_chains'))$metric
bimodal_trajectory_df <- 
  rbind(bimodal_trajectory_df, 
        data.frame(t=c(1:chain_length), type='indep_coupling', trajectory=5, 
                   metric=indep_coupling, metric_sd=NA))

# save(bimodal_trajectory_df, file="stylized_example_ula_mala/bimodal_trajectory_df.RData")
# load("stylized_example_ula_mala/bimodal_trajectory_df.RData")

plot.labels <- unname(TeX(c('$CUB_{1,t}$ (I=1)','$CUB_{1,t}$ (I=1000)','True $W_1$', 'Indep. Coupling')))
bimodal_trajectory_plot1a <- 
  ggplot(bimodal_trajectory_df %>% filter(type != '1_chain_ergodic', t <= 1000), 
         aes(x=t, y=metric, color=interaction(type, trajectory), 
             linetype=interaction(type, trajectory))) + 
  geom_line(size=1) +
  geom_ribbon(aes(ymin=metric-2*metric_sd, ymax=metric+2*metric_sd, 
                  fill=as.factor(trajectory)),alpha=0.2, colour = NA) +
  xlab(TeX('iteration $t$')) + ylab(TeX('$W_1$ upper bounds')) +
  scale_colour_manual(name=TeX(''), labels=plot.labels,
                      values = c('lightgrey','black','black', 'black','black')) +
  scale_linetype_manual(name=TeX(''), labels=plot.labels,
                        values = c('dotted','solid', 'dashed', 'dotted', 'dotdash')) +
  scale_fill_manual(name=TeX(''), labels=plot.labels,
                    values = c('NA','black', 'NA', 'NA', 'NA')) +
  #scale_linetype_discrete(guide = "none") +
  #scale_linetype_manual(breaks = c("coupling", "coupling_averaged", "true_w1", "indep_coupling"))+
  scale_x_continuous(breaks=seq(0,1e3,5e2)) +
  scale_y_continuous(breaks=seq(0,12,2), limits = c(0,6)) +
  theme_classic(base_size = 16) + 
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
bimodal_trajectory_plot1a
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/implementation/bimodal_trajectory_plot1a.pdf", plot = bimodal_trajectory_plot1a, width = 4, height = 5)

plot.labels <- unname(TeX(c('$CUB_{1,t}$ (I=1)','$CUB_{1}$ (I=1)', 'True $W_1$', 'Indep. Coupling')))
bimodal_trajectory_plot1b <- 
  ggplot(bimodal_trajectory_df %>% filter(type != '100_chains', t <= Inf), 
         aes(x=t, y=metric, color=interaction(type, trajectory), 
             linetype=interaction(type, trajectory))) + 
  geom_line(size=1) +
  geom_ribbon(aes(ymin=metric-2*metric_sd,ymax=metric+2*metric_sd, fill=as.factor(trajectory)),
              alpha=0.2, colour = NA) +
  xlab(TeX('iteration $t$')) + ylab(TeX('$W_1$ upper bounds')) +
  scale_colour_manual(name=TeX(''), labels=plot.labels,
                      values = c('lightgrey','black','black', 'black','black')) +
  scale_linetype_manual(name=TeX(''), labels=plot.labels,
                        values = c('dotted','solid', 'dashed', 'dotted', 'dotdash')) +
  scale_fill_manual(name=TeX(''), labels=plot.labels,
                    values = c('NA','black', 'NA', 'NA', 'NA')) +
  #scale_linetype_discrete(guide = "none") +
  #scale_linetype_manual(breaks = c("coupling", "coupling_averaged", "true_w1", "indep_coupling"))+
  scale_x_continuous(breaks=seq(0,1e4,1e3)) +
  scale_y_continuous(breaks=seq(0,12,2), limits = c(0,6)) +
  theme_classic(base_size = 16) + 
  theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
bimodal_trajectory_plot1b
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/implementation/bimodal_trajectory_plot1b.pdf", plot = bimodal_trajectory_plot1b, width = 4, height = 5)







