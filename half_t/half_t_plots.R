# Half-t approximate MCMC plots

# load("half_t/half_t_trajectory_df.RData")

library(dplyr)
library(ggplot2)
library(latex2exp)

# setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
# load("half_t/half_t_bounds_df_riboflavin.RData")
# load("half_t/half_t_bounds_df_synthetic.RData")

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

