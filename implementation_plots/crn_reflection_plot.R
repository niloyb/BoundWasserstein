rm(list = ls())
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source(file = 'kernels.R')
source(file = 'estimators.R')
source(file = 'implementation_plots/mvn_mixture_functions.R')

# Metric considered: L1
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}

################################################################################
#  Targets: mixture of gaussians
dimension <- 1
no_components_P <- 2
weights_P <- rep(1,no_components_P)
weights_P <- weights_P/sum(weights_P)
means_P <- cbind(rep(-2,dimension),rep(2,dimension))

weights_Q <- rep(1,no_components_P)
weights_Q <- weights_Q/sum(weights_Q)
means_Q <- cbind(rep(-1,dimension),rep(1,dimension))

# # # Marginal MALA chains
chain_length <- 5000
init_state_P <- SamplerP(1)[1,]
sigma_mh_P <- 2
mala_output_P <- mala(init_state_P, LogPdfP, GradLogPdfP, sigma_mh_P, chain_length)
mean(mala_output_P[2:chain_length,1]!=mala_output_P[1:(chain_length-1),1])
init_state_Q <- SamplerQ(1)[1,]
sigma_mh_Q <- 2
mala_output_Q <- mala(init_state_Q, LogPdfQ, GradLogPdfQ, sigma_mh_Q, chain_length)
mean(mala_output_Q[2:chain_length,1]!=mala_output_Q[1:(chain_length-1),1])

par(mfrow=c(3,2))
matplot(mala_output_P[,1], type = 'l')
matplot(mala_output_Q[,1], type = 'l')
hist(mala_output_P[,1], 100)
hist(mala_output_Q[,1], 100)
hist(SamplerP(chain_length)[,1], 100)
hist(SamplerQ(chain_length)[,1], 100)

# CRN Coupled MALA
chain_length <- 1000
sigma_mh <- c(sigma_mh_P,sigma_mh_Q)

coupled_crn_mala_sampler <- 
  function(){
    return(coupled_mala(SamplerP(1)[1,], SamplerQ(1)[1,], LogPdfP, LogPdfQ,
                        GradLogPdfP, GradLogPdfQ, sigma_mh, chain_length))
  }
coupled_reflect_mala_sampler <- 
  function(){
    return(coupled_mala(SamplerP(1)[1,], SamplerQ(1)[1,], LogPdfP, LogPdfQ,
                        GradLogPdfP, GradLogPdfQ, sigma_mh, chain_length,
                        reflect_threshold=Inf))
  }

p <- 1
crn_ula_mala_trajectory <- 
  wp_ub_estimate(coupled_crn_mala_sampler, no_chains=1000,
                 p=p, metric=metric_l2, parallel=TRUE)
reflect_ula_mala_trajectory <- 
  wp_ub_estimate(coupled_reflect_mala_sampler, no_chains=1000, p=p, 
                 metric=metric_l2, parallel=TRUE)
crn_reflection_df <-
  rbind(data.frame(t=c(1:chain_length), type='crn', trajectory=1,
                   metric=(crn_ula_mala_trajectory$wp_power_p_ub_mean_tracjectory)^(1/p)),
        data.frame(t=c(1:chain_length), type='reflect', trajectory=2,
                   metric=(reflect_ula_mala_trajectory$wp_power_p_ub_mean_tracjectory)^(1/p)))

crn_reflection_bounds_df <-
  rbind(data.frame(type='crn', ub=crn_ula_mala_trajectory$wp_ub, ub_se=crn_ula_mala_trajectory$wp_ub_se),
        data.frame(type='reflect', ub=reflect_ula_mala_trajectory$wp_ub, ub_se=reflect_ula_mala_trajectory$wp_ub_se))


# True W1L1 distance
# If dimension = 1
Q_samples <- SamplerQ(1000)
P_samples <- SamplerP(1000)
true_w1 <- mean(((sort(Q_samples)-sort(P_samples))^2)^0.5)
# For dimension >= 2
# Q_samples <- SamplerQ(1000)
# Q_samples <- transport::pp(Q_samples)
# P_samples <- SamplerP(1000)
# P_samples <- transport::pp(P_samples)
# true_w1 <- transport::wasserstein(Q_samples, P_samples)
crn_reflection_df <- 
  rbind(crn_reflection_df, 
        data.frame(t=c(1:chain_length), type='true_w1', trajectory=3, metric=true_w1))

crn_reflection_bounds_df <- 
  rbind(crn_reflection_bounds_df, 
        data.frame(type='true_w1', ub=true_w1, ub_se=0))

# Independent coupling
# indep_coupling <- mean(metric_l2(SamplerQ(1000), SamplerP(1000)))
indep_coupling_mean <- mean(metric_l2(SamplerQ(1000), SamplerP(1000))) # When dimension=1
indep_coupling_sd <- sd(metric_l2(SamplerQ(1000), SamplerP(1000)))/sqrt(1000)
# indep_coupling <- # As we initialize the chain from the indep coupling
#   (crn_reflection_df %>% filter(t==1, type == 'crn'))$metric
crn_reflection_df <- 
  rbind(crn_reflection_df, 
        data.frame(t=c(1:chain_length), type='indep_coupling', trajectory=4, metric=indep_coupling_mean))

crn_reflection_bounds_df <- 
  rbind(crn_reflection_bounds_df, 
        data.frame(type='indep_coupling', ub=indep_coupling_mean, ub_se=indep_coupling_sd))

plot.labels <- unname(TeX(c('$CUB_{1,t}$ (CRN)','$CUB_{1,t}$ (Reflect)','True $W_1$', 'Indep. Coupling')))
# plot.labels <- unname(TeX(c('CRN Coupling','Reflection Coupling', '$W_1$ Optimal', 'Independent')))
crn_reflection_plot <- 
  ggplot(crn_reflection_df, 
         aes(x=t, y=metric, color=interaction(type, trajectory), 
             linetype=interaction(type, trajectory))) + 
  geom_line(size=1) +
  xlab(TeX('iteration $t$')) + ylab(TeX('$W_1$ upper bounds')) +
  scale_colour_manual(name=TeX(''), labels=plot.labels,
                      values = c('lightgrey', rep('black',3))) +
  scale_linetype_manual(name=TeX(''), labels=plot.labels,
                        values = c('solid', 'solid', 'dotted','dotdash')) +
  #scale_linetype_discrete(guide = "none") +
  #scale_linetype_manual(breaks = c("coupling", "coupling_averaged", "true_w1", "indep_coupling"))+
  # scale_x_continuous(breaks=seq(0,1e3,100), limits = c(0,200)) +
  scale_y_continuous(breaks=seq(0,12,1), limits = c(0,2.5)) +
  theme_classic(base_size = 16) + 
  theme(legend.position = 'bottom', legend.key.width=unit(0.9,"cm")) +
  guides(color=guide_legend(nrow=2,byrow=TRUE))
crn_reflection_plot
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/implementation/crn_reflection_plot.pdf", plot = crn_reflection_plot, width = 4, height = 5)

