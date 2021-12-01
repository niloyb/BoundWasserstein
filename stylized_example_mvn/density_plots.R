# Bivariate density plot multivariate normals

# Libraries
library(mvtnorm)
library(dplyr)
library(ggplot2)
library(latex2exp)
library(gridExtra)


# Set parameters of Gaussians
d <- 2
mu1 <- rep(0, d)
mu2 <- rep(0, d)
sigma1 <- diag(d)
rho <- 0.5
sigma2 <- matrix(c(1,rho,rho,1),2,2)

# Generate pdf
x <- seq(-3, 3, by=0.01)
y <- seq(-3, 3, by=0.01)

pdf1_df <- data.frame(expand.grid(x, y))
colnames(pdf1_df) <- c('x','y')
pdf2_df <- pdf1_df
pdf1_df <- pdf1_df %>% 
  dplyr::group_by(x,y) %>%
  dplyr::mutate(pdf = mvtnorm::dmvnorm(c(x,y), mean = mu1, sigma = sigma1),
                distribution='P')
pdf2_df <- pdf2_df %>% 
  dplyr::group_by(x,y) %>%
  dplyr::mutate(pdf = mvtnorm::dmvnorm(c(x,y), mean = mu2, sigma = sigma2),
                distribution='Q')
pdf_df <- rbind(pdf1_df,pdf2_df)

# Plots
# plot1 <- ggplot(pdf1_df, aes(x,y,z=pdf)) + 
#   stat_contour(geom="polygon", aes(fill=-..level..)) + 
#   # stat_density2d(aes(fill = ..density..), geom = "raster", contour = FALSE)+
#   # geom_raster(aes(fill = -pdf)) +
#   theme_classic(base_size = 18) + 
#   xlab(TeX('$x$')) + ylab(TeX('$y$')) +
#   scale_fill_gradientn(colours = gray.colors(2)) +
#   theme(legend.position = "none") + coord_cartesian(xlim=c(-3,3), ylim=c(-3,3))
# 
# plot2 <- ggplot(pdf2_df, aes(x, y, z = pdf)) +
#   stat_contour(geom="polygon", aes(fill=-..level..)) + 
#   # geom_raster(aes(fill = -pdf)) +
#   theme_classic(base_size = 18) + 
#   xlab(TeX('$x$')) + ylab(TeX('$y$')) +
#   scale_fill_gradientn(colours = gray.colors(2)) +
#   theme(legend.position = "none") +
#   coord_cartesian(xlim=c(-3,3),
#                   ylim=c(-3,3))

plot3 <- ggplot(pdf_df %>% dplyr::group_by(distribution), 
                aes(x, y, z = pdf)) +
  stat_contour(geom="polygon", aes(fill=distribution), alpha=0.5) + 
  # xlab(TeX('$x$')) + ylab(TeX('$y$')) +
  xlab(element_blank()) +ylab(element_blank()) +
  scale_fill_manual(name="Distribution",
                    breaks = c("P", "Q"),
                    labels = unname(TeX(c("$N(0,I)$", "$N(0,\\Sigma)$"))),
                    values = grey.colors(2),
                    guide = "legend") +
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3)) +
  theme_classic(base_size = 16)
plot3
# ggsave(filename = "/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/writeup/images/stylized_example1/bivariate_plot.pdf", plot = plot3, width = 4, height = 3)
