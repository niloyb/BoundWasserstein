# Half-t approximate MCMC plots

# load("half_t/half_t_trajectory_df.RData")

library(dplyr)
library(ggplot2)
library(latex2exp)

# trajectories_df_plot <- trajectories_df %>% 
#   dplyr::group_by(approxMCMCparam, t_dist_df, t) %>%
#   dplyr::summarise(metric=mean(metric))
# 
# W2L2_bounds_df <- trajectories_df %>% 
#   dplyr::group_by(approxMCMCparam, t_dist_df) %>%
#   dplyr::summarise(W2L2UB=mean(W2L2UB), W2L2LB=mean(W2L2LB))
# 
# 
# #
# half_t_trajectory_plot <- 
#   ggplot(trajectories_df_plot %>% filter(approxMCMCparam < 0.1, t>0), 
#          aes(x=t, y=metric, color=as.factor(approxMCMCparam),
#              linetype=as.factor(t_dist_df))) +
#   # aes(x=t, y=metric, color=as.factor(model), linetype=as.factor(model))) + 
#   geom_line(size=1) + xlab(TeX('iteration')) + 
#   ylab(TeX('$W_2^2$ upper bound')) +
#   theme_classic(base_size = 14) + 
#   scale_colour_manual(name=TeX('Threshold $\\epsilon$'), 
#                       values = c('black', 'grey')) +
#   scale_y_continuous(trans = 'log10') +
#   theme(legend.position = 'bottom', legend.key.width=unit(1,"cm")) +
#   guides(color=guide_legend(nrow=2))+
#   guides(linetype=guide_legend(nrow=2, title=TeX('$\\nu$')))
# half_t_trajectory_plot
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/half_t/half_t_trajetory.pdf", plot = half_t_trajectory_plot, width = 4, height = 4)
# 
# wass_bounds_plot <- 
#   ggplot(W2L2_bounds_df %>% filter(approxMCMCparam < 0.1), 
#          aes(x=as.factor(approxMCMCparam), y=W2L2UB)) + 
#   # geom_point(size=4, aes(shape=as.factor(t_dist_df))) + 
#   # xlab(TeX('')) + 
#   geom_errorbar(aes(ymin=W2L2LB, ymax=W2L2UB, linetype=as.factor(t_dist_df)),
#                 width=.5, position=position_dodge(.9)) +
#   ylab(TeX('$W_2$ upper and lower bounds')) +
#   xlab(TeX('Threshold $\\epsilon$')) +
#   #scale_x_continuous(trans = 'log10') +
#   #scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
#   #scale_y_continuous(limits=c(0,3)) + 
#   theme_classic(base_size = 14) +
#   theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
#   guides(linetype=guide_legend(TeX('$\\nu$'))) + theme(legend.position = 'bottom')
# wass_bounds_plot
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/half_t/wass_bounds_plot.pdf", plot = wass_bounds_plot, width = 4, height = 4)
# 


### NEW PLOTS
# setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
# load("half_t/half_t_bounds_df_maize.RData")
# load("half_t/half_t_bounds_df_synthetic.RData")

# half_t_wass_bounds_plot <- 
#   ggplot(bounds_df, aes(x=as.factor(approxMCMCparam), y=W2L2UBmean)) + 
#   geom_crossbar(aes(ymin=W2L2UBmean-W2L2UBsd, 
#                     ymax=W2L2UBmean+W2L2UBsd,
#                     group=as.factor(t_dist_df)), 
#                 # width=10, 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_crossbar(aes(ymin=W2L2LBmean-W2L2LBsd, 
#                     ymax=W2L2LBmean+W2L2LBsd,
#                     group=as.factor(t_dist_df)), 
#                 # width=10, 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_errorbar(aes(ymin=W2L2LBmean, ymax=W2L2UBmean, linetype=as.factor(t_dist_df)), 
#                 #width=.5, 
#                 position=position_dodge(.9)) +
#   # geom_crossbar(aes(ymin=indep_W2L2UBmean-indep_W2L2UBsd, 
#   #                   ymax=indep_W2L2UBmean+indep_W2L2UBsd,
#   #                   linetype=as.factor(t_dist_df)), 
#   #               # width=10, 
#   #               color='grey', fill='grey',
#   #               position=position_dodge(.9)) +
#   # geom_point(size=7, aes(y=indep_W2L2UBmean), shape=4, color=as.factor(t_dist_df)) + 
#   ylab(TeX('$W_2$ upper and lower bounds')) +
#   xlab(TeX('Approx. MCMC threshold $\\epsilon$')) +
#   # scale_x_continuous(labels = scales::scientific) +
#   theme_classic(base_size = 14) +
#   guides(linetype=guide_legend(TeX('Half-t prior degree of freedom $\\nu$'))) + 
#   theme(legend.position = 'bottom')
# half_t_wass_bounds_plot
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/half_t/half_t_wass_bounds_plot.pdf", plot = half_t_wass_bounds_plot, width = 6, height = 4)



half_t_wass_bounds_plot <- 
  ggplot(bounds_df, aes(x=as.factor(approxMCMCparam), y=W2L2UBmean)) + 
  geom_crossbar(aes(ymin=W2L2UBmean-W2L2UBsd, 
                    ymax=W2L2UBmean+W2L2UBsd,
                    group=as.factor(t_dist_df)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean-W2L2LBsd, 
                    ymax=W2L2LBmean+W2L2LBsd,
                    group=as.factor(t_dist_df)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2LBmean, ymax=W2L2UBmean, group=as.factor(t_dist_df)), 
                #width=.5, 
                position=position_dodge(.9)) +
  # geom_crossbar(aes(ymin=indep_W2L2UBmean-indep_W2L2UBsd, 
  #                   ymax=indep_W2L2UBmean+indep_W2L2UBsd,
  #                   linetype=as.factor(t_dist_df)), 
  #               # width=10, 
  #               color='grey', fill='grey',
  #               position=position_dodge(.9)) +
  # geom_point(size=7, aes(y=indep_W2L2UBmean), shape=4, color=as.factor(t_dist_df)) + 
  ylab(TeX('$W_2$ upper and lower bounds')) +
  xlab(TeX('Approx. MCMC threshold $\\epsilon$')) +
  # scale_x_continuous(labels = scales::scientific) +
  # scale_y_continuous(limits = c(0,0.8)) +
  theme_classic(base_size = 14) +
  guides(linetype=guide_legend(TeX('Half-t prior degree of freedom $\\nu$'))) + 
  theme(legend.position = 'bottom')
half_t_wass_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/half_t/half_t_wass_bounds_riboflavin_plot.pdf", plot = half_t_wass_bounds_plot, width = 4, height = 3)
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/half_t/half_t_wass_bounds_synthetic_plot.pdf", plot = half_t_wass_bounds_plot, width = 4, height = 3)

