# MVN mixture functions

library(MASS)
library(mvtnorm)

## MVN mixture
# weights_P, means_P globally defined
# weights_P is a vector of length no. of components
# means_P is a matrix of size dim * no. of components
# LogPdfP <- function(x){
#   weights_P <- weights_P/sum(weights_P)
#   return(log(sum(weights_P*exp(colSums(-0.5*(means_P-x)^2)))))
# }
LogPdfP <- function(x){
  weights_P <- weights_P/sum(weights_P)
  M <- max(colSums(-0.5*(means_P-x)^2))
  normalized_vector <- colSums(-0.5*(means_P-x)^2)-M
  return(M + log(sum(weights_P*exp(normalized_vector))))
}
# Gradient of log density
GradLogPdfP <- function(x){
  weights_P <- weights_P/sum(weights_P)
  # weighted unnormalized pdf by component
  M <- max(colSums(-0.5*(means_P-x)^2))
  normalized_vector <- colSums(-0.5*(means_P-x)^2)-M
  weighted_pdf_by_component <- weights_P*exp(normalized_vector)
  weighted_pdf_by_component <- 
    weighted_pdf_by_component/sum(weighted_pdf_by_component)
  return(rowSums(t(t(means_P-x)*weighted_pdf_by_component)))
}
SamplerP <-function(no_indep_samples){
  weights_P <- weights_P/sum(weights_P)
  sampled_compoments <- 
    sample(length(weights_P), no_indep_samples, 
           prob = weights_P, replace = TRUE)
  output <- 
    matrix(rnorm(dimension*no_indep_samples), 
           nrow = no_indep_samples, ncol = dimension)
  for (i in c(1:no_indep_samples)){
    output[i,] <- output[i,] + means_P[,sampled_compoments[i]]
  }
  return(output)
}

# weights_Q, means_Q globally defined
# weights_Q is a vector of length no. of components
# means_Q is a matrix of size dim * no. of components
# LogPdfQ <- function(x){
#   weights_Q <- weights_Q/sum(weights_Q)
#   return(log(sum(weights_Q*exp(colSums(-0.5*(means_Q-x)^2)))))
# }
LogPdfQ <- function(x){
  weights_Q <- weights_Q/sum(weights_Q)
  M <- max(colSums(-0.5*(means_Q-x)^2))
  normalized_vector <- colSums(-0.5*(means_Q-x)^2)-M
  return(M + log(sum(weights_Q*exp(normalized_vector))))
}
# Gradient of log density
GradLogPdfQ <- function(x){
  weights_Q <- weights_Q/sum(weights_Q)
  # weighted unnormalized pdf by component
  M <- max(colSums(-0.5*(means_Q-x)^2))
  normalized_vector <- colSums(-0.5*(means_Q-x)^2)-M
  weighted_pdf_by_component <- weights_Q*exp(normalized_vector)
  weighted_pdf_by_component <- 
    weighted_pdf_by_component/sum(weighted_pdf_by_component)
  return(rowSums(t(t(means_Q-x)*weighted_pdf_by_component)))
}
SamplerQ <-function(no_indep_samples){
  weights_Q <- weights_Q/sum(weights_Q)
  sampled_compoments <- 
    sample(length(weights_Q), no_indep_samples, 
           prob = weights_Q, replace = TRUE)
  output <- 
    matrix(rnorm(dimension*no_indep_samples), 
           nrow = no_indep_samples, ncol = dimension)
  for (i in c(1:no_indep_samples)){
    output[i,] <- output[i,] + means_Q[,sampled_compoments[i]]
  }
  return(output)
}



