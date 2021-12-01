rm(list = ls())
source('../huggins_comparison/viabel_functions.R')
source('../paper_examples/stylized_example_mvn/mvn_functions.R')
source('../paper_examples/rcpp_functions.R')

# Test 1: estimator_is function 
# remember IS is abysmally useless if target far away from proposal ...
integrand_function <- function(x){return(x)}

weight_function <- function(x){
  return(exp(dnorm(x, mean=2, log = TRUE)-dnorm(x, mean=0, log = TRUE)))
}
samples <- rnorm(1000000)
estimates_is <- estimator_is(integrand_function, weight_function, samples)
mean(estimates_is)
sd(estimates_is)

# Test 2a: CUBO and ELBO functions
log_M <- 2
unnormalised_target_logpdf <- function(x){dnorm(x,mean=1,log=TRUE)+log_M}
proposal_logpdf <- function(x){dnorm(x,log=TRUE)}
proposal_samples <- t(as.matrix(rnorm(100000), nrow=1))

# Should get estimator_cubo >= log_M >= estimator_elbo
estimator_cubo(unnormalised_target_logpdf, proposal_logpdf, proposal_samples)
estimator_elbo(unnormalised_target_logpdf, proposal_logpdf, proposal_samples)



# Test 2b: CUBO and ELBO functions with MVN case
dimension <- 100
alpha <- 0.5
# Covariance matrix
cov_mat <- cov_matrix(alpha, dimension)
inverse_cov_matrix <- cov_matrix_inv(alpha, dimension)

LogPdfP_normalized <- 
  function(x) {
    # dimension <- length(x)
    norm_constant <- -(dimension/2)*log(2*pi)-(1/2)*log(det(cov_mat))
    return(norm_constant+LogPdfP(x))
  }

LogPdfQ_normalized <- 
  function(x) {
    # dimension <- length(x)
    norm_constant <- -(dimension/2)*log(2*pi)
    return(norm_constant+LogPdfQ(x))
  }

unnormalised_target_logpdf <- LogPdfP
# NOTE: need proposal_logpdf to be normalized
proposal_logpdf <- LogPdfQ_normalized 
proposal_samples <- SamplerQ(10000)

# Should get estimator_cubo >= log_M >= estimator_elbo
estimator_cubo(unnormalised_target_logpdf, proposal_logpdf, proposal_samples)
estimator_elbo(unnormalised_target_logpdf, proposal_logpdf, proposal_samples)
# Here log_M 
log_M <- LogPdfP(proposal_samples[,1,drop=FALSE])-LogPdfP_normalized(proposal_samples[,1,drop=FALSE])
log_M






# Test 3: W2L2 upper bound
no_samples <- 100000
proposal_samples <- SamplerQ(no_samples)
proposal_componentwise_variances <- diag(cov_mat)

elbo_proposal_logpdf <- LogPdfP_normalized
elbo_proposal_samples <- SamplerP(no_samples)

# elbo_proposal_logpdf <- LogPdfQ
# elbo_proposal_samples <- SamplerQ(no_samples)


# Should equal log_M as elbo_proposal_logpdf=LogPdfP_normalized
# estimator_elbo(unnormalised_target_logpdf, elbo_proposal_logpdf, elbo_proposal_samples)
# log_M
# 0.5*(dimension*log(2*pi)+log(det(cov_mat)))

# Should be bigger than log_M
# estimator_cubo(unnormalised_target_logpdf, proposal_logpdf, proposal_samples)


cubo_elbo_diff <- 
  estimator_cubo(unnormalised_target_logpdf, proposal_logpdf, proposal_samples) - 
  estimator_elbo(unnormalised_target_logpdf, elbo_proposal_logpdf, elbo_proposal_samples)

estimator_h(unnormalised_target_logpdf, cubo_proposal_logpdf=proposal_logpdf, 
            elbo_proposal_logpdf=LogPdfP_normalized, 
            cubo_proposal_samples=proposal_samples, elbo_proposal_samples=elbo_proposal_samples)


huggins_W2L2_ub_gaussian(unnormalised_target_logpdf, proposal_logpdf, proposal_samples, 
           proposal_componentwise_variances, elbo_proposal_logpdf, elbo_proposal_samples)
# Independent coupling based upper bound
(2*dimension)^0.5

cov_mat_sqrt <- expm::sqrtm(cov_mat)
true_W2_squared_mean <- norm(cov_mat_sqrt-diag(dimension), type = 'F')^2
true_W2_squared_mean^0.5




componentwise_variances <- proposal_componentwise_variances



