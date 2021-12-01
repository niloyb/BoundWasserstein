# Functions for the Bayesian Logistic Regression Example
# Note: response takes values in {-1,1}

# Stable Logsigmoid function
logSigmoid <- function(x){
  # pos_x <- x>0
  # x[pos_x] <- -log(1+exp(-x[pos_x]))
  # x[!pos_x] <- x[!pos_x]-log(1+exp(x[!pos_x]))
  # return(x)
  return(-log(1+exp(-abs(x)))+(x-abs(x))/2)
}

# Log target pdf up to a constant
# X, y, prior_mean, inv_prior_cov_matrix globally defined
targetLogPdf <- function(beta){
  log_prior <- -0.5*car::wcrossprod((beta-prior_mean), (beta-prior_mean), inv_prior_cov_matrix)
  log_likelihood <- sum(logSigmoid((y*(X%*%beta))))
  return(log_prior + log_likelihood)
}

# Grad of Log target pdf
# X, y, prior_mean, inv_prior_cov_matrix globally defined
targetGradLogPdf <- function(beta){
  grad_log_likelihood <- colSums(X*((exp(logSigmoid(-y*(X%*%beta)))*y)[,1]))
  grad_log_prior <- -(crossprod(inv_prior_cov_matrix, (beta-prior_mean)))[,1]
  return(grad_log_prior+grad_log_likelihood)
}

## Approx distribtion N(mu, Sigma) 
# Log approx distribution pdf up to a constant
# mu, inv_sigma globally defined
mvnLogPdf <- function(beta) {
  return(-0.5*car::wcrossprod((beta-mu), (beta-mu), inv_sigma))
}

# Grad of Log proposal pdf
mvnGradLogPdf <- function(beta) {
  return(-crossprod(inv_sigma, (beta-mu))[,1])
}

# Log target pdf up to a constant
# X, y, prior_mean, inv_prior_cov_matrix, subsample_ratio globally defined
approxTargetLogPdf <- function(beta){
  n <- dim(X)[1]
  sampled_data <- sample(c(1:n), (n*subsample_ratio), replace = FALSE)
  X <- X[sampled_data,]
  y <- y[sampled_data]
  log_prior <- -0.5*car::wcrossprod((beta-prior_mean), (beta-prior_mean), inv_prior_cov_matrix)
  log_likelihood <- sum(logSigmoid((y*(X%*%beta))))/subsample_ratio
  return(log_prior + log_likelihood)
}

# Grad of Log target pdf
# X, y, prior_mean, inv_prior_cov_matrix, subsample_ratio globally defined
approxTargetGradLogPdf <- function(beta){
  n <- dim(X)[1]
  sampled_data <- sample(c(1:n), (n*subsample_ratio), replace = FALSE)
  X <- X[sampled_data,]
  y <- y[sampled_data]
  grad_log_likelihood <- colSums(X*((exp(logSigmoid(-y*(X%*%beta)))*y)[,1]))/subsample_ratio
  grad_log_prior <- -(crossprod(inv_prior_cov_matrix, (beta-prior_mean)))[,1]
  return(grad_log_prior+grad_log_likelihood)
}


