# Testing functions for the Bayesian Logistic Regression Example

# Logistic Regression Wasserstein Distance Plots
rm(list = ls())
set.seed(1)
source(file='/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/bayesian_logistic_regression/logreg_functions.R')
source(file='/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/kernels.R')
source(file='/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/rcpp_functions.R')


# Generating data
n <- 100000
p <- 10
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
s <- 5
true_beta <- matrix(0,p,1)
true_beta[1:s] = 2^(-(seq(s)/4-9/4))
probs <- exp(logSigmoid(X%*%true_beta))
y <- 2*rbinom(n, size = 1, prob = probs)-1
# Initialising and defining first chain: posterior of logistic regression
beta1_init <- rnorm(p)
prior_mean <- rep(0,p)
inv_prior_cov_matrix <- diag(p)
sigma_mh <- 0.1

# Starting beta1_init from (approximately) stationarity
beta1_long_run <- ula(beta1_init, targetGradLogPdf, sigma_mh, 100)

beta1_long_run2 <- mala(beta1_init, targetLogPdf, targetGradLogPdf, sigma_mh, 100)


beta1_init <- beta1_long_run[10000,]

# Initialising and defining second chain: N(mu, sigma)
#mu <- true_beta
#inv_sigma <- diag(p)
mu <- colMeans(beta1_long_run[c(5e3:1e4),])
inv_sigma <- diag(1/apply(beta1_long_run[c(5e3:1e4),], 2, var))
# Starting beta1_init from stationarity
beta2_init <- mvtnorm::rmvnorm(1, mean = mu, sigma = solve(inv_sigma))[1,]

# Running one chain (with synchronous coupling)
chain_length <- 1000
coupled_chain <- coupled_ula(beta2_init, beta1_init, mvnGradLogPdf, targetGradLogPdf,
                             sigma_mh, chain_length)
# Marginal traceplots
par(mfrow=c(2,1))
matplot(coupled_chain$q[,c(1:10)], type = 'l')
matplot(coupled_chain$pi[,c(1:10)], type = 'l')
# Plotting distance
par(mfrow=c(1,1))
matplot(rowSums((coupled_chain$q-coupled_chain$pi)^2)^0.5, type = 'l')
par(mfrow=c(1,1))
matplot(abs(coupled_chain$q-coupled_chain$pi)[,c(15)], type = 'l')

