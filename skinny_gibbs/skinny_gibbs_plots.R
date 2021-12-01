############## Generating Plots for the HMC example ###########

library(dplyr)
library(ggplot2)
library(latex2exp)


#setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
load("skinny_gibbs/skinny_gibbs_bounds_df.RData")
wass_beta_bounds_plot <- 
  ggplot(skinny_gibbs_df, aes(x=dimension, y=W2L2UBmean_beta)) + 
  geom_crossbar(aes(ymin=W2L2UBmean_beta-W2L2UBsd_beta/sqrt(no_chains), 
                    ymax=W2L2UBmean_beta+W2L2UBsd_beta/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean_beta-W2L2LBsd_beta/sqrt(no_chains), 
                    ymax=W2L2LBmean_beta+W2L2LBsd_beta/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2UBmean_beta, ymax=W2L2LBmean_beta), 
                #width=.5, 
                position=position_dodge(.9)) +
  ylab(TeX('$W_2$ upper and lower bounds on $\\beta$')) + xlab(TeX('Synthetic dataset dimension $d$')) + 
  scale_y_continuous(limits=c(0,3)) +
  theme_classic(base_size = 14) # +
#theme(legend.position = 'bottom') + guides(point=guide_legend("Parameterization", nrow=2,byrow=TRUE))
wass_beta_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_bounds_plot_beta.pdf", plot = wass_beta_bounds_plot, width = 4, height = 3.5)

wass_z_bounds_plot <-
  ggplot(skinny_gibbs_df, aes(x=dimension, y=W2L2UBmean_z)) +
  geom_crossbar(aes(ymin=W2L2UBmean_z-W2L2UBsd_z/sqrt(no_chains),
                    ymax=W2L2UBmean_z+W2L2UBsd_z/sqrt(no_chains)),
                # width=10,
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean_z-W2L2LBsd_z/sqrt(no_chains),
                    ymax=W2L2LBmean_z+W2L2LBsd_z/sqrt(no_chains)),
                # width=10,
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2UBmean_z, ymax=W2L2LBmean_z),
                #width=.5,
                position=position_dodge(.9)) +
  ylab(TeX('$W_2$ upper and lower bounds on $Z$')) + xlab(TeX('Synthetic dataset dimension $d$')) +
  scale_y_continuous(limits=c(0,4)) +
  theme_classic(base_size = 14) # +
#theme(legend.position = 'bottom') + guides(point=guide_legend("Parameterization", nrow=2,byrow=TRUE))
wass_z_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_bounds_plot_z.pdf", plot = wass_z_bounds_plot, width = 4, height = 3.5)


################################################################################
load("skinny_gibbs/skinny_gibbs_genome_df.RData")

skinny_gibbs_genome_df$dataset[1] <- "1.malware"
skinny_gibbs_genome_df$dataset[2] <- "2.lymph"
plot.breaks <- c("1.malware", "2.lymph")
plot.labels <- c('Malware', 'Lymph Node')

wass_beta_bounds_genome_plot <- 
  ggplot(skinny_gibbs_genome_df, aes(x=dataset, y=W2L2UBmean_beta)) + 
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  geom_crossbar(aes(ymin=W2L2UBmean_beta-W2L2UBsd_beta/sqrt(no_chains), 
                    ymax=W2L2UBmean_beta+W2L2UBsd_beta/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean_beta-W2L2LBsd_beta/sqrt(no_chains), 
                    ymax=W2L2LBmean_beta+W2L2LBsd_beta/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2UBmean_beta, ymax=W2L2LBmean_beta), 
                #width=.5, 
                position=position_dodge(.9)) +
  ylab(TeX('$W_2$ upper and lower bounds')) + xlab(TeX('Datasets')) + 
  scale_y_continuous(limits=c(0,10)) +
  theme_classic(base_size = 18) # +
# theme(legend.position = 'bottom') + guides(point=guide_legend("Parameterization", nrow=2,byrow=TRUE))
wass_beta_bounds_genome_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_beta_bounds_genome_plot.pdf", plot = wass_beta_bounds_genome_plot, width = 8, height = 4)


wass_z_bounds_genome_plot <- 
  ggplot(skinny_gibbs_genome_df, aes(x=dataset, y=W2L2UBmean_z)) + 
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  geom_crossbar(aes(ymin=W2L2UBmean_z-W2L2UBsd_z/sqrt(no_chains), 
                    ymax=W2L2UBmean_z+W2L2UBsd_z/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean_z-W2L2LBsd_z/sqrt(no_chains), 
                    ymax=W2L2LBmean_z+W2L2LBsd_z/sqrt(no_chains)), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_errorbar(aes(ymin=W2L2UBmean_z, ymax=W2L2LBmean_z), 
                #width=.5, 
                position=position_dodge(.9)) +
  ylab(TeX('$W_2$ upper and lower bounds on $Z$')) + xlab(TeX('Datasets')) + 
  scale_y_continuous(limits=c(0,8)) +
  theme_classic(base_size = 20) # +
# theme(legend.position = 'bottom') + guides(point=guide_legend("Parameterization", nrow=2,byrow=TRUE))
wass_z_bounds_genome_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_z_bounds_genome_plot.pdf", plot = wass_z_bounds_genome_plot, width = 8, height = 4)





# #setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
# load("skinny_gibbs/skinny_gibbs_df.RData")
# wass_bounds_plot <- 
#   ggplot(skinny_gibbs_df, aes(x=p, y=coupled_W1hammingUBmean)) + 
#   geom_crossbar(aes(ymin=coupled_W1hammingUBmean-coupled_W1hammingUBsd/sqrt(no_chains), 
#                     ymax=coupled_W1hammingUBmean+coupled_W1hammingUBsd/sqrt(no_chains)), 
#                 # width=10, 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_crossbar(aes(ymin=W1hammingLBmean-W1hammingLBsd/sqrt(no_chains), 
#                     ymax=W1hammingLBmean+W1hammingLBsd/sqrt(no_chains)), 
#                 # width=10, 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_crossbar(aes(ymin=indep_W1hammingUBmean-indep_W1hammingUBsd/sqrt(no_chains), 
#                     ymax=indep_W1hammingUBmean+indep_W1hammingUBsd/sqrt(no_chains)), 
#                 # width=10, 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_errorbar(aes(ymin=W1hammingLBmean, ymax=coupled_W1hammingUBmean), 
#                 #width=.5, 
#                 position=position_dodge(.9)) +
#   # geom_hline(aes(yintercept=indep_W1hammingUBmean, linetype='dotted')) +
#   geom_point(size=7, aes(y=indep_W1hammingUBmean), shape=4) + 
#   ylab(TeX('$W_1$ upper and lower bounds')) + xlab(TeX('Dimension $d$')) + 
#   scale_y_continuous(limits=c(0,20)) +
#   theme_classic(base_size = 14) # +
#   #theme(legend.position = 'bottom') + guides(point=guide_legend("Parameterization", nrow=2,byrow=TRUE))
# wass_bounds_plot
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_bounds_plot_z.pdf", plot = wass_bounds_plot, width = 4, height = 3.5)
# 
# 
# 
# wass_bounds_plot_gwas <- 
#   ggplot(skinny_gibbs_genome_df, aes(x=dataset, y=coupled_W1hammingUBmean)) + 
#   geom_crossbar(aes(ymin=coupled_W1hammingUBmean-coupled_W1hammingUBsd/sqrt(no_chains), 
#                     ymax=coupled_W1hammingUBmean+coupled_W1hammingUBsd/sqrt(no_chains)), 
#                 # width=0.5,
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_crossbar(aes(ymin=W1hammingLBmean-W1hammingLBsd/sqrt(no_chains), 
#                     ymax=W1hammingLBmean+W1hammingLBsd/sqrt(no_chains)), 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_crossbar(aes(ymin=indep_W1hammingUBmean-indep_W1hammingUBsd/sqrt(no_chains), 
#                     ymax=indep_W1hammingUBmean+indep_W1hammingUBsd/sqrt(no_chains)), 
#                 color='grey', fill='grey',
#                 position=position_dodge(.9)) +
#   geom_errorbar(aes(ymin=W1hammingLBmean, ymax=coupled_W1hammingUBmean), 
#                 #width=.5, 
#                 position=position_dodge(.9)) +
#   geom_point(size=7, aes(y=indep_W1hammingUBmean), shape=4) + 
#   ylab(TeX('$W_1$ upper and lower bounds')) + xlab(TeX('Datasets')) + 
#   scale_y_continuous(limits=c(0,80)) +
#   theme_classic(base_size = 14) # +
# #theme(legend.position = 'bottom') + guides(point=guide_legend("Parameterization", nrow=2,byrow=TRUE))
# wass_bounds_plot_gwas
# # ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_bounds_plot_z_gwas.pdf", plot = wass_bounds_plot_gwas, width = 4, height = 3.5)



