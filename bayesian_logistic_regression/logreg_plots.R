############## Generating Plots for the logistic regression example ###########
rm(list = ls())
set.seed(runif(1))

library(ggplot2)
library(dplyr)
library(latex2exp)

#setwd('Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
load("bayesian_logistic_regression/logreg_bounds_df_pima.RData")
logreg_df_pima <- bounds_df
load("bayesian_logistic_regression/logreg_bounds_df2_ds1.RData")
logreg_df_ds1 <- bounds_df

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
  ylab(TeX('$W_2$ upper and lower bounds')) + xlab(TeX('')) + 
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  scale_y_continuous(limits=c(0,0.1))
ds1_bounds_plot
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/ds1_bounds_plot.pdf", plot = ds1_bounds_plot, width = 4, height = 4)


# Combined plot
logreg_df_combined <- rbind(logreg_df_pima, logreg_df_ds1)
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
  xlab(TeX('Approximate MCMC or variational procedure')) + 
  ylab(TeX('$W_2$ upper and lower bounds')) +
  scale_x_discrete(breaks=plot.breaks, labels=plot.labels) +
  # scale_y_continuous(trans = 'log10') +
  scale_y_continuous(limits =c(0,0.3)) +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +
  guides(linetype=guide_legend("Dataset")) + theme(legend.position = 'right')
combined_ub_lb_plot
#ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/logistic_regression/combined_ub_lb_plot.pdf", plot = combined_ub_lb_plot, width = 6, height = 4)


