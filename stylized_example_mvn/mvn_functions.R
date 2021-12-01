# Couplings to obtain Wasserstein distance upper bounds between two Gaussians

library(MASS)
library(mvtnorm)

# Sample from Q (dimension globally defined)
SamplerQ <- function(no_indep_samples){
  return(t(mvtnorm::rmvnorm(n = no_indep_samples, mean = rep(0,dimension))))
}
# Q log density up to normalizing constant
LogPdfQ <- function(x) {
  log_pdf <- -0.5*sum(x^2)
  return(log_pdf)
}
# Gradient of log density
GradLogPdfQ <- function(x) {
  return(-x)
}
## Distribution P
# P is mean zero MVN with decaying non-diag cov matrix
# Covariance matrix of P
cov_matrix <- function(alpha, dimension){
  Sigma_pi <- 0.5*diag(1, dimension, dimension)
  if(dimension>1){
    for (i in 2:dimension){
      for (j in 1:(i-1)){
        Sigma_pi[i,j] <- alpha^(abs(i-j))
      }
    }
  }
  Sigma_pi <- Sigma_pi + t(Sigma_pi)
  return(Sigma_pi)
}
cov_matrix_inv <- function(alpha, dimension){
  return(chol2inv(chol(cov_matrix(alpha, dimension))))
}
# Sample from P (dimension globally defined)
SamplerP <-function(no_indep_samples){
  return(t(mvtnorm::rmvnorm(n = no_indep_samples, mean = rep(0,dimension), sigma = cov_mat)))
}
# P log density up to normalizing constant
LogPdfP <- function(x) {
  #return(-0.5*sum((x%*%inverse_cov_matrix)*x))
  return(-0.5*sum(x*cpp_prod(inverse_cov_matrix, matrix(x, ncol = 1))))
}
# Gradient of log density
GradLogPdfP <- function(x){
  # return(-inverse_cov_matrix%*%x)
  return(-cpp_prod(inverse_cov_matrix, matrix(x, ncol = 1)))
}




