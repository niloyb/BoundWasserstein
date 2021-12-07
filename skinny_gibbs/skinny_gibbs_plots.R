############## Generating Plots for the HMC example ###########

library(dplyr)
library(ggplot2)
library(latex2exp)


# setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
# load("skinny_gibbs/skinny_gibbs_df.RData")

skinny_gibbs_df$dataset[1] <- "1.malware"
skinny_gibbs_df$dataset[2] <- "2.lymph"
plot.breaks <- c("1.malware", "2.lymph")
plot.labels <- c('Malware', 'Lymph Node')

wass_beta_bounds_plot <- 
  ggplot(skinny_gibbs_df, aes(x=dataset, y=W2L2UBmean_beta)) + 
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  geom_crossbar(aes(ymin=W2L2UBmean_beta-W2L2UBsd_beta, 
                    ymax=W2L2UBmean_beta+W2L2UBsd_beta), 
                # width=10, 
                color='grey', fill='grey',
                position=position_dodge(.9)) +
  geom_crossbar(aes(ymin=W2L2LBmean_beta-W2L2LBsd_beta, 
                    ymax=W2L2LBmean_beta+W2L2LBsd_beta), 
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
wass_beta_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/skinny_gibbs/wass_beta_bounds_plot.pdf", plot = wass_beta_bounds_plot, width = 8, height = 4)

