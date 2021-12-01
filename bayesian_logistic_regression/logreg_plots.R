############## Generating Plots for the logistic regression example ###########
rm(list = ls())
set.seed(runif(1))

library(ggplot2)
library(dplyr)
library(latex2exp)

#setwd('Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
#load("bayesian_logistic_regression/logreg_trajectory_df_pima.RData")
load("bayesian_logistic_regression/logreg_bounds_df_pima.RData")
logreg_df_pima <- bounds_df
load("bayesian_logistic_regression/logreg_bounds_df_ds1.RData")
logreg_df_ds1 <- bounds_df

# # Plot 1: Pima trajectory
# # Wasserstein distance trajectories between MALA and AVDI, Laplace Approximations
# plot.breaks <- c("mf_advi", "laplace", "sgld0.1", "sgld0.5", "ula")
# plot.labels <- unname(TeX(c("Mean Field VB", "Laplace", 'SGLD10%','SGLD50%','ULA')))
# plot.name <- TeX('')
# logreg_plot <- 
#   ggplot(logreg_df_pima$trajectory %>% filter(t <=300), 
#          aes(x=t, y=metric, color=as.factor(type), linetype=as.factor(type))) + 
#   geom_line(size=1) + xlab(TeX('iteration')) + 
#   ylab(TeX('$W_2^2$ upper bound')) +
#   scale_colour_manual(name=plot.name, breaks = plot.breaks, labels= plot.labels,
#                       values = c('darkgrey', 'darkgrey', 'black','black', 'black')) +
#   scale_linetype_manual(name=plot.name, breaks = plot.breaks, labels=plot.labels,
#                         values = c('dashed', 'solid', 'dotted','dashed', 'solid')) +
#   theme_classic(base_size = 14) + 
#   scale_y_continuous(trans = 'log10') + 
#   # scale_y_continuous(breaks=seq(0,10,0.1), limits=c(0,0.5)) + 
#   theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) + 
#   guides(color=guide_legend(nrow=3,byrow=TRUE))
# logreg_plot
# #ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/logreg_plot_pima.pdf", plot = logreg_plot, width = 4, height = 5)

# # Plot 2: ds1 trajectory
# plot.breaks <- c("mf_advi", "laplace", "sgld0.1", "sgld0.5", "ula")
# plot.labels <- unname(TeX(c("Mean Field VB", "Laplace", 'SGLD10%','SGLD50%','ULA')))
# plot.name <- TeX('')
# logreg_plot <- 
#   ggplot(logreg_df_ds1$trajectory %>% filter(t <=300), 
#          aes(x=t, y=metric, color=as.factor(type), linetype=as.factor(type))) + 
#   geom_line(size=1) + xlab(TeX('iteration')) + 
#   ylab(TeX('$W_2^2$ upper bound')) +
#   scale_colour_manual(name=plot.name, breaks = plot.breaks, labels= plot.labels,
#                       values = c('darkgrey', 'darkgrey', 'black','black', 'black')) +
#   scale_linetype_manual(name=plot.name, breaks = plot.breaks, labels=plot.labels,
#                         values = c('dashed', 'solid', 'dotted','dashed', 'solid')) +
#   theme_classic(base_size = 14) + 
#   scale_y_continuous(breaks=10^seq(-8,-2,2), trans = 'log10') + 
#   # scale_y_continuous(breaks=seq(0,10,0.1), limits=c(0,0.5)) + 
#   theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) + 
#   guides(color=guide_legend(nrow=3,byrow=TRUE))
# logreg_plot
# #ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/logreg_plot_ds1.pdf", plot = logreg_plot, width = 4, height = 5)

# ordered.labels.pima <- paste(c(2,5,1,3,4),logreg_df_pima$marginal_errors$approx, sep = '')
# logreg_df_pima$marginal_errors$approx <- ordered.labels.pima
# logreg_df_pima$marginal_errors$dataset <- 'Pima'
# ordered.labels.ds1 <- paste(c(2,5,1,3,4),logreg_df_ds1$marginal_errors$approx, sep = '')
# logreg_df_ds1$marginal_errors$approx <- ordered.labels.ds1
# logreg_df_ds1$marginal_errors$dataset <- 'DS1'
# marginal_errors_df <- rbind(logreg_df_pima$marginal_errors, logreg_df_ds1$marginal_errors)
# 
# # Plot 3: Marginal Covariance Error
# plot.breaks <- c("2mf_advi", "5laplace", "1sgld0.1", "3sgld0.5", "4ula")
# plot.labels <- unname(TeX(c("Mean Field VB", "Laplace", 'SGLD 10%','SGLD 50%','ULA')))
# plot.name <- TeX('Approx.')
# cov_error_plot <- 
#   ggplot(marginal_errors_df, aes(x=as.factor(approx), y=relaive_cov_error)) + 
#   geom_point(size=4, aes(shape=as.factor(dataset))) + xlab(TeX('')) + 
#   ylab(TeX('Covariance Error $| \ \\Sigma^{-1} \\hat{\\Sigma} - I_d \ |_2$')) +
#   scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
#   scale_y_continuous(limits=c(0.1,4), trans = 'log10') +
#   theme_classic(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   guides(shape=guide_legend("Dataset")) + theme(legend.position = 'bottom')
# cov_error_plot
# #ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/cov_error_plot.pdf", plot = cov_error_plot, width = 4, height = 5)

# wass_bounds_plot <- 
#   ggplot(marginal_errors_df, aes(x=as.factor(approx), y=W2L2UB)) + 
#   # geom_point(size=4, aes(shape=as.factor(dataset))) + 
#   xlab(TeX('')) + 
#   geom_errorbar(aes(ymin=W2L2LB, ymax=W2L2UB, linetype=as.factor(dataset)),
#                 width=.5, position=position_dodge(.9)) +
#   ylab(TeX('$W_2$ upper and lower bounds')) +
#   scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
#   scale_y_continuous(limits=c(0,0.3)) +
#   theme_classic(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   guides(linetype=guide_legend("Dataset")) + theme(legend.position = 'bottom')
# wass_bounds_plot
# #ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/wass_bounds_plot.pdf", plot = wass_bounds_plot, width = 4, height = 5)








# no_chains <- 100
# Pima Wassertein Bounds plot
ordered.labels.pima <- paste(c(2,5,1,3,4),logreg_df_pima$type, sep = '')
logreg_df_pima$type <- ordered.labels.pima
logreg_df_pima$dataset <- 'Pima'

plot.breaks <- c("2mf_advi", "5laplace", "1sgld0.1", "3sgld0.5", "4ula")
plot.labels <- unname(TeX(c("Mean Field VB", "Laplace", 'SGLD 10%','SGLD 50%','ULA')))
plot.name <- TeX('Approx.')
pima_bounds_plot <- 
  ggplot(logreg_df_pima, aes(x=as.factor(type), y=W2L2UBmean)) + 
  geom_crossbar(aes(ymin=W2L2UBmean-W2L2UBsd/sqrt(no_chains), 
                    ymax=W2L2UBmean+W2L2UBsd/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean-W2L2LBsd/sqrt(no_chains), 
                    ymax=W2L2LBmean+W2L2LBsd/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2LBmean, ymax=W2L2UBmean), 
                #width=.5, 
                position=position_dodge(.9)) +
  # geom_hline(aes(yintercept=indep_W1hammingUBmean, linetype='dotted')) +
  # geom_crossbar(aes(ymin=indep_W2L2UBmean-indep_W2L2UBsd/sqrt(no_chains), 
  #                   ymax=indep_W2L2UBmean+indep_W2L2UBsd/sqrt(no_chains)), 
  #               # width=10, 
  #               color='grey', fill='grey',
  #               position=position_dodge(.9)) +
  # geom_point(size=7, aes(y=indep_W2L2UBmean), shape=4) + 
  ylab(TeX('$W_2$ upper and lower bounds')) + xlab(TeX('')) + 
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_y_continuous(limits=c(0,0.5))
pima_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/pima_bounds_plot.pdf", plot = pima_bounds_plot, width = 4, height = 4)
  
# DS1 Wassertein Bounds plot
ordered.labels.pima <- paste(c(2,5,1,3,4),logreg_df_ds1$type, sep = '')
logreg_df_ds1$type <- ordered.labels.pima
logreg_df_ds1$dataset <- 'DS1'

plot.breaks <- c("2mf_advi", "5laplace", "1sgld0.1", "3sgld0.5", "4ula")
plot.labels <- unname(TeX(c("Mean Field VB", "Laplace", 'SGLD 10%','SGLD 50%','ULA')))
plot.name <- TeX('Approx.')
ds1_bounds_plot <- 
  ggplot(logreg_df_ds1, aes(x=as.factor(type), y=W2L2UBmean)) + 
  geom_crossbar(aes(ymin=W2L2UBmean-W2L2UBsd/sqrt(no_chains), 
                    ymax=W2L2UBmean+W2L2UBsd/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean-W2L2LBsd/sqrt(no_chains), 
                    ymax=W2L2LBmean+W2L2LBsd/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2LBmean, ymax=W2L2UBmean), 
                #width=.5, 
                position=position_dodge(.9)) +
  # geom_hline(aes(yintercept=indep_W1hammingUBmean, linetype='dotted')) +
  # geom_crossbar(aes(ymin=indep_W2L2UBmean-indep_W2L2UBsd/sqrt(no_chains), 
  #                   ymax=indep_W2L2UBmean+indep_W2L2UBsd/sqrt(no_chains)), 
  #               # width=10, 
  #               color='grey', fill='grey',
  #               position=position_dodge(.9)) +
  # geom_point(size=7, aes(y=indep_W2L2UBmean), shape=4) + 
  ylab(TeX('$W_2$ upper and lower bounds')) + xlab(TeX('')) + 
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_y_continuous(limits=c(0,0.1))
ds1_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/ds1_bounds_plot.pdf", plot = ds1_bounds_plot, width = 4, height = 4)


# Combined plot
logreg_df_combined <- rbind(logreg_df_pima, logreg_df_ds1)

# Plot 3: Together
plot.breaks <- c("2mf_advi", "5laplace", "1sgld0.1", "3sgld0.5", "4ula")
plot.labels <- unname(TeX(c("Mean Field VB", "Laplace", 'SGLD 10%','SGLD 50%','ULA')))
plot.name <- TeX('Approx.')
combined_ub_plot <- 
  ggplot(logreg_df_combined, aes(x=as.factor(type), y=W2L2UBmean)) + 
  geom_point(size=4, aes(shape=as.factor(dataset))) + 
  xlab(TeX('Sampling Algorithm')) + 
  ylab(TeX('$W_2$ upper and lower bounds')) +
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  scale_y_continuous(trans = 'log10') +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  guides(shape=guide_legend("Dataset")) + theme(legend.position = 'bottom')
combined_ub_plot
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/combined_ub_plot.pdf", plot = combined_ub_plot, width = 4, height = 5)



# ggplot(logreg_df_combined, aes(x=as.factor(type), y=W2L2UBmean)) + 
#   geom_point(size=4, aes(shape=as.factor(dataset)), position=position_dodge(.9))

combined_ub_lb_plot <- 
  ggplot(logreg_df_combined, aes(x=as.factor(type), y=W2L2UBmean)) + 
  geom_crossbar(aes(x=as.factor(type), ymin=W2L2UBmean-W2L2UBsd/sqrt(no_chains), 
                    ymax=W2L2UBmean-W2L2UBsd/sqrt(no_chains),group=as.factor(dataset)), 
                fill='grey', color='gray', position=position_dodge(.9)) +
  geom_crossbar(aes(x=as.factor(type), ymin=W2L2LBmean-W2L2LBsd/sqrt(no_chains), 
                    ymax=W2L2LBmean-W2L2LBsd/sqrt(no_chains),group=as.factor(dataset)), 
                fill='grey', color='gray', position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2LBmean, ymax=W2L2UBmean, linetype=as.factor(dataset)), 
                position=position_dodge(.9)) +
  xlab(TeX('Sampling Algorithm')) + 
  ylab(TeX('$W_2$ upper and lower bounds')) +
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  # scale_y_continuous(trans = 'log10') +
  scale_y_continuous(limits =c(0,0.3)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  guides(linetype=guide_legend("Dataset")) + theme(legend.position = 'right')
combined_ub_lb_plot
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/combined_ub_lb_plot.pdf", plot = combined_ub_lb_plot, width = 6, height = 4)


