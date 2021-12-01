# Functions to implement the Wasserstein distance upper bounds of
# Validated Variational Inference via Practical Posterior Error Bounds. 
# Huggins, Kasprzak, Campbell and Broderick. AISTATS, 2020.

estimator_is <- function(integrand_function, weight_function, samples, normalize=FALSE){
  weights <- weight_function(samples)
  if(normalize){weights <- weights/sum(weights)}
  return(integrand_function(samples)*weights)
}

estimator_cubo <- function(unnormalised_target_logpdf, proposal_logpdf, proposal_samples){
  weight_function <- function(x){exp(unnormalised_target_logpdf(x)-proposal_logpdf(x))}
  no_samples <- dim(proposal_samples)[2]
  weights <- apply(proposal_samples, 2, weight_function)
  return(0.5*log(mean(weights^2)))
}

estimator_elbo <- function(unnormalised_target_logpdf, proposal_logpdf, proposal_samples){
  log_weight_function <- function(x){unnormalised_target_logpdf(x)-proposal_logpdf(x)}
  no_samples <- dim(proposal_samples)[2]
  log_weights <- apply(proposal_samples, 2, log_weight_function)
  return(mean(log_weights))
}

estimator_h <- 
  function(unnormalised_target_logpdf, cubo_proposal_logpdf, elbo_proposal_logpdf, 
           cubo_proposal_samples, elbo_proposal_samples){
    cubo_elbo_diff <- 
      estimator_cubo(unnormalised_target_logpdf, cubo_proposal_logpdf, cubo_proposal_samples) - 
      estimator_elbo(unnormalised_target_logpdf, elbo_proposal_logpdf, elbo_proposal_samples)
    return(2*max(cubo_elbo_diff,0))
}

constant_ei_ind_gaussian <- function(componentwise_variances){
  # dimension <- length(componentwise_variances)
  ub_function <- function(epsilon){
    return((3/2-(1/2)*sum(log(1-2*epsilon*componentwise_variances)))/epsilon)
  }
  grid <- seq(0.0001, 1/(2*max(componentwise_variances))-0.0001, 0.0001)
  values <- sapply(grid, ub_function)
  return(2*min(values)^0.5)
}

# Update to allow componentwise_variances with non-unique values
constant_pi_ind_gaussian <- function(componentwise_variances){
  return(2*(2*sum(componentwise_variances^2)+sum(componentwise_variances)^2)^(1/4))
}


huggins_W2L2_ub_gaussian <- 
  function(unnormalised_target_logpdf, proposal_logpdf, proposal_samples, 
           proposal_componentwise_variances, elbo_proposal_logpdf, elbo_proposal_samples){
    constant_ei <- constant_ei_ind_gaussian(proposal_componentwise_variances)
    estimate_h <- 
      estimator_h(unnormalised_target_logpdf, cubo_proposal_logpdf=proposal_logpdf, 
                  elbo_proposal_logpdf=elbo_proposal_logpdf, 
                  cubo_proposal_samples=proposal_samples,
                  elbo_proposal_samples=elbo_proposal_samples)
    output_ei <- constant_ei*(estimate_h^(1/2)+(estimate_h/2)^(1/4))
    constant_pi <- constant_pi_ind_gaussian(proposal_componentwise_variances)
    output_pi <- constant_pi*(exp(estimate_h)-1)^(1/4)
    return(min(output_ei,output_pi,na.rm=TRUE))
    # return(output_ei)
  }










