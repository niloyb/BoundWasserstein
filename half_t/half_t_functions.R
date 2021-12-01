################ Bayesian High dimensional regression functions ################
# Code modified from the code corresponding to the article
# "Coupled Markov chain Monte Carlo for high-dimensional regression with Half-t priors" 
# (https://arxiv.org/abs/2012.04798)

################ Single chain MCMC functions ################
# Functions for a blocked Gibbs sampler with half-t priors
require(Rcpp)
# require(inline)
require(RcppEigen)

## Rcpp functions ##
Rcpp::cppFunction('
Eigen::MatrixXd cpp_crossprod(const Eigen::MatrixXd X){
Eigen::MatrixXd output;
output = (X.transpose())*X;
return output;
}', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
Eigen::MatrixXd output;
output = X*Y;
return output;
}', depends = 'RcppEigen')

## Matrix calculation helper functions ##
X_eta_tX <- function(eta, X, X_transpose){
  # NOTE: (1) For MacOS with veclib BLAS, crossprod is fast via multi-threading
  return(crossprod(X_transpose*c(1/eta)^0.5))
  # return(cpp_crossprod(X_transpose*c(1/eta)^0.5))
}
M_matrix <- function(xi, eta, X_eta_tX_matrix, n){
  if(length(eta)==0) return(diag(n))
  return(diag(n) + (xi^-1)*X_eta_tX_matrix)
}

## xi update given eta ##
# Unnormalized posterior pdf of log(xi)
log_ratio <- function(xi, eta, X_eta_tX_matrix, y, a0, b0)
{
  n <- length(y)
  M <- M_matrix(xi,eta,X_eta_tX_matrix,n)
  chol_M <- chol(M)
  log_det_M <- 2*sum(log(diag(chol_M)))
  M_inverse <- chol2inv(chol_M)
  ssr <- b0 + t(y)%*%((M_inverse)%*%y)
  log_likelihood <- -0.5*log_det_M -0.5*(n+a0)*log(ssr)
  log_prob <- -log(sqrt(xi)*(1+xi))
  return(list('log_likelihood'=log_likelihood+log_prob,
              'ssr' = ssr, 'M_matrix_inverse' = M_inverse))
}
# Unnormalized posterior pdf of log(xi) for approximate MCMC
log_ratio_approx <- function(xi, eta, X, X_transpose, y, a0, b0, active_set)
{
  n <- length(y)
  log_prob <- -log(sqrt(xi)*(1+xi))
  if (sum(active_set)==0)
  {
    M_inverse <- diag(n)
    ssr <- b0 + sum(y^2)
    log_likelihood <- -0.5*(n+a0)*log(ssr)
  } else{
    eta <- eta[active_set]
    X <- X[,active_set,drop=F]
    X_transpose <- X_transpose[active_set,,drop=F]
    if(sum(active_set)==1){
      woodbury_matrix_part <- xi*eta + cpp_prod(X_transpose, X)
    } else{
      woodbury_matrix_part <- xi*diag(eta) + cpp_prod(X_transpose, X)
    }
    woodbury_matrix_part_inverse <- chol2inv(chol(woodbury_matrix_part))
    M_inverse <- diag(n) -
      cpp_prod(cpp_prod(X,woodbury_matrix_part_inverse),X_transpose)
    log_det_M <- sum( log( xi^(-1)*(svd((eta^(-0.5))*X_transpose)$d)^2 + 1 ) )
    ssr <- b0 + t(y)%*%((M_inverse)%*%y)
    log_likelihood <- -0.5*log_det_M -0.5*(n+a0)*log(ssr)
  }
  return(list('log_likelihood'=log_likelihood+log_prob, 'ssr' = ssr,
              'M_matrix_inverse' = M_inverse))
}
# Metropolis-Hastings update of xi given eta
xi_update <- function(current_xi, eta, X, X_transpose, y, a0, b0, std_MH,
                      approximate_algo_delta=0, fixed=FALSE)
{
  n <- length(y)
  if (fixed==TRUE){
    min_xi <- current_xi
    active_set <- ((min_xi*eta)^(-1) > approximate_algo_delta)
    
    # Matrix calculations
    if (sum(active_set)>n) # Woodbury inversion when sum(active_set)>n
    {
      X_eta_tX_matrix <- X_eta_tX(eta[active_set], X[ ,active_set,drop=F],
                                  X_transpose[active_set, ,drop=F])
      log_ratio_current_and_ssr <- log_ratio(current_xi, eta[active_set],
                                             X_eta_tX_matrix, y, a0, b0)
    } else { # Woodbury inversion when sum(active_set)<n
      log_ratio_current_and_ssr <-
        log_ratio_approx(current_xi, eta, X, X_transpose, y, a0, b0, active_set)
    }
  } else {
    proposed_xi <- exp(rnorm(1, log(current_xi), std_MH))
    min_xi <- min(current_xi, proposed_xi)
    active_set <- ((min_xi*eta)^(-1) > approximate_algo_delta)
    
    # Matrix calculations
    if (sum(active_set)>n) # Woodbury inversion when sum(active_set)>n
    {
      X_eta_tX_matrix <- X_eta_tX(eta[active_set],X[ ,active_set,drop=F],
                                  X_transpose[active_set, ,drop=F])
      log_ratio_current_and_ssr <- log_ratio(current_xi, eta[active_set],
                                             X_eta_tX_matrix, y, a0, b0)
      log_ratio_proposed_and_ssr <- log_ratio(proposed_xi, eta[active_set],
                                              X_eta_tX_matrix, y, a0, b0)
    } else { # Woodbury inversion when sum(active_set)<n
      log_ratio_current_and_ssr <-
        log_ratio_approx(current_xi, eta, X, X_transpose, y, a0, b0, active_set)
      log_ratio_proposed_and_ssr <-
        log_ratio_approx(proposed_xi, eta, X, X_transpose, y, a0, b0, active_set)
    }
    
    # MH accept- reject step
    log_accept_prob <-
      (log_ratio_proposed_and_ssr$log_likelihood -
         log_ratio_current_and_ssr$log_likelihood) +
      (log(proposed_xi)-log(current_xi))
    
    if (log(runif(1))<log_accept_prob)
    {
      current_xi <- proposed_xi
      log_ratio_current_and_ssr <- log_ratio_proposed_and_ssr
    }
  }
  return(list('xi'=current_xi, 'ssr' = log_ratio_current_and_ssr$ssr,
              'M_matrix' = log_ratio_current_and_ssr$M_matrix_inverse,
              'active_set' = active_set))
}

## sigma^2 update step given eta and xi
sigma2_update <- function(xi, eta, n, ssr, a0, b0, sigma2_fixed_value=NULL)
{
  if(!is.null(sigma2_fixed_value)) return(sigma2_fixed_value)
  return(1/(rgamma(1, shape = (n+a0)/2, rate = (ssr)/2)))
}

## beta update given xi, eta, sigma2
beta_update <- function(xi, sigma2, eta, X, X_transpose, y,
                        M_matrix_inverse, active_set)
{
  p <- length(eta)
  n <- length(y)
  u = rnorm(p, 0, 1)
  u = (sqrt(xi*eta)^-1)*u
  v = X%*%u + c(rnorm(n,0,1))
  v_star <- M_matrix_inverse%*%(y/sqrt(sigma2) - v)
  U <- (xi^-1)* ( ((eta[active_set])^(-1))*(X_transpose[active_set, ,drop=F]) )
  u[active_set] <- u[active_set] + U%*%v_star
  beta <- sqrt(sigma2)*(u)
  return(beta)
}


## Full blocked Gibbs samplers ##
half_t_kernel <-
  function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           xi_current, sigma2_current,
           beta_current, eta_current, approximate_algo_delta,
           nrepeats_eta = 1, verbose = FALSE,
           xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df)
  {
    if (verbose) ptm <- proc.time()
    n <- dim(X)[1]
    p <- dim(X)[2]
    # xi update
    xi_new <- xi_update(xi_current, eta_current, X, X_transpose, y, a0, b0, std_MH,
                        approximate_algo_delta, fixed=xi_fixed)
    # sigma2 update
    sigma2_fixed_value <- NULL
    if(sigma2_fixed==TRUE) sigma2_fixed_value <- sigma2_current
    sigma2_new <- sigma2_update(xi_new$xi, eta_current, n, xi_new$ssr, a0, b0,
                                sigma2_fixed_value)
    # beta update
    beta_new <- beta_update(xi_new$xi, sigma2_new, eta_current, X, X_transpose, y,
                            xi_new$M_matrix, xi_new$active_set)
    # eta update
    eta_new <- eta_update_half_t(xi_new$xi, sigma2_new, beta_new, eta_current,
                                 t_dist_df, nrepeats_eta)
    
    if (verbose) print(proc.time()[3]-ptm[3])
    output <- list( 'beta_samples'=beta_new, 'eta_samples'=eta_new,
                    'sigma2_samples'=sigma2_new, 'xi_samples'=xi_new$xi)
    return(output)
  }

#' half_t_mcmc
half_t_mcmc <-
  function(chain_length, burnin, X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           rinit=NULL, approximate_algo_delta=0, nrepeats_eta = 1,
           verbose = FALSE, xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df)
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    if(is.null(rinit)){
      # Initializing from the prior
      rinit <- function(){
        xi <- (1/rt(1, df=1))^2
        sigma2 <- 1/rgamma(1, shape = a0/2, rate = b0/2)
        eta <- (1/rt(p, df=t_dist_df))^2
        beta <- rnorm(p)*sqrt(sigma2/(xi*eta))
        return(list(xi = xi, sigma2 = sigma2, beta = beta, eta = eta))
      }
    }
    # Initializing first chain
    chain <- rinit()
    xi_samples <- rep(NA, (chain_length-burnin))
    sigma2_samples <- rep(NA, (chain_length-burnin))
    beta_samples <- matrix(NA, nrow=(chain_length-burnin), ncol=p)
    eta_samples <- matrix(NA, nrow=(chain_length-burnin), ncol=p)
    xi_samples[1] <- xi_current <- chain$xi
    sigma2_samples[1] <- sigma2_current <- chain$sigma2
    beta_samples[1,] <- beta_current <- chain$beta
    eta_samples[1,] <- eta_current <- chain$eta
    
    i <- 1
    while (i <= chain_length)
    {
      output <-
        half_t_kernel(X, X_transpose, y, a0=a0, b0=b0, std_MH=std_MH,
                      xi_current, sigma2_current, beta_current, eta_current,
                      approximate_algo_delta, nrepeats_eta = nrepeats_eta,
                      verbose = verbose, xi_fixed=xi_fixed,
                      sigma2_fixed=sigma2_fixed, t_dist_df)
      xi_current <- output$xi_samples
      sigma2_current <- output$sigma2_samples
      beta_current <- output$beta_samples
      eta_current <- output$eta_samples
      if(i > burnin)
      {
        xi_samples[(i-burnin)] <- xi_current
        sigma2_samples[(i-burnin)] <- sigma2_current
        beta_samples[(i-burnin),] <- beta_current
        eta_samples[(i-burnin),] <- eta_current
      }
      if (verbose) print(i)
      i <- i + 1
    }
    return(list( 'beta_samples'=beta_samples, 'eta_samples'=eta_samples,
                 'sigma2_samples'=sigma2_samples, 'xi_samples'=xi_samples))
  }




################ Coupled chain MCMC functions ################
# Functions for the coupled blocked Gibbs sampler with half-t priors
require(lubridate)

## Couplings for xi update given eta ##
# Coupled Metropolis-Hastings update of xi given eta
crn_max_xi_coupling <- 
  function(current_xi_1, eta_1, current_xi_2, eta_2,
           X, X_transpose, y, a0, b0, std_MH,
           approximate_algo_delta_1=0, approximate_algo_delta_2=0,
           epsilon_xi=0, fixed=FALSE)
{
  if(fixed==TRUE){ # When xi is fixed in the Gibbs sampler
    min_xi_1 <- current_xi_1
    min_xi_2 <- current_xi_2
    active_set_1 <- ((min_xi_1*eta_1)^(-1) > approximate_algo_delta_1)
    active_set_2 <- ((min_xi_2*eta_2)^(-1) > approximate_algo_delta_2)
    
    if (sum(active_set_1)>n)
    {
      X_eta_tX_matrix_1 <-
        X_eta_tX(eta_1[active_set_1],X[, active_set_1, drop=F],
                 X_transpose[active_set_1, , drop=F])
      log_ratio_current_ssr_matrixinv_1 <-
        log_ratio(current_xi_1, eta_1[active_set_1], X_eta_tX_matrix_1, y, a0, b0)
    } else {
      log_ratio_current_ssr_matrixinv_1 <-
        log_ratio_approx(current_xi_1, eta_1, X, X_transpose, y, a0, b0, active_set_1)
    }
    
    if (sum(active_set_2)>n)
    {
      X_eta_tX_matrix_2 <-
        X_eta_tX(eta_2[active_set_2],X[, active_set_2, , drop=F], X_transpose[active_set_2, , drop=F])
      log_ratio_current_ssr_matrixinv_2 <-
        log_ratio(current_xi_2, eta_2[active_set_2], X_eta_tX_matrix_2, y, a0, b0)
    } else {
      log_ratio_current_ssr_matrixinv_2 <-
        log_ratio_approx(current_xi_2, eta_2, X, X_transpose, y, a0, b0, active_set_2)
    }
  } else { # When xi is varying in the Gibbs sampler
    standard_normal <- rnorm(1, mean = 0, sd = 1)
    log_proposed_xi_1 <- standard_normal*sqrt(std_MH) + log(current_xi_1)
    
    relative_error_delta <- abs(log(current_xi_1)-log(current_xi_2))
    # print(relative_error_delta)
    if ((0 < relative_error_delta) & (relative_error_delta < epsilon_xi)) {
      # using max coupling to get the proposal in the MH-algo
      if( dnorm(log_proposed_xi_1, mean = log(current_xi_1), sd = std_MH, log = TRUE) +
          log(runif(1)) < dnorm(log_proposed_xi_1, mean = log(current_xi_2), sd = std_MH, log = TRUE) )
      {
        log_proposed_xi_2 <- log_proposed_xi_1
      } else {
        reject <- TRUE
        y_proposal <- NA
        attempts <- 0
        while(reject){
          attempts <- attempts + 1
          # print(attempts)
          y_proposal <- rnorm(1, mean = log(current_xi_2), sd = std_MH)
          reject <- ( dnorm(y_proposal, mean = log(current_xi_2), sd = std_MH, log = TRUE) +
                        log(runif(1)) < dnorm(y_proposal, mean = log(current_xi_1), sd = std_MH, log = TRUE) )
        }
        log_proposed_xi_2 <- y_proposal
      }
    } else {
      # using common random numbers to get the proposal in the MH-algo
      log_proposed_xi_2 <- standard_normal*sqrt(std_MH) + log(current_xi_2)
    }
    
    proposed_xi_1 <- exp(log_proposed_xi_1)
    proposed_xi_2 <- exp(log_proposed_xi_2)
    
    min_xi_1 <- min(current_xi_1, proposed_xi_1)
    min_xi_2 <- min(current_xi_2, proposed_xi_2)
    active_set_1 <- ((min_xi_1*eta_1)^(-1) > approximate_algo_delta_1)
    active_set_2 <- ((min_xi_2*eta_2)^(-1) > approximate_algo_delta_2)
    
    if (sum(active_set_1)>n)
    {
      X_eta_tX_matrix_1 <- X_eta_tX(eta_1[active_set_1],X[, active_set_1, drop=F], X_transpose[active_set_1, , drop=F])
      log_ratio_current_ssr_matrixinv_1 <- log_ratio(current_xi_1, eta_1[active_set_1], X_eta_tX_matrix_1, y, a0, b0)
      log_ratio_proposed_ssr_matrixinv_1 <- log_ratio(proposed_xi_1, eta_1[active_set_1], X_eta_tX_matrix_1, y, a0, b0)
    } else {
      log_ratio_current_ssr_matrixinv_1 <- log_ratio_approx(current_xi_1, eta_1, X, X_transpose, y, a0, b0, active_set_1)
      log_ratio_proposed_ssr_matrixinv_1 <- log_ratio_approx(proposed_xi_1, eta_1, X, X_transpose, y, a0, b0, active_set_1)
    }
    
    if (sum(active_set_2)>n)
    {
      X_eta_tX_matrix_2 <- X_eta_tX(eta_2[active_set_2],X[, active_set_2, , drop=F], X_transpose[active_set_2, , drop=F])
      log_ratio_current_ssr_matrixinv_2 <- log_ratio(current_xi_2, eta_2[active_set_2], X_eta_tX_matrix_2, y, a0, b0)
      log_ratio_proposed_ssr_matrixinv_2 <- log_ratio(proposed_xi_2, eta_2[active_set_2], X_eta_tX_matrix_2, y, a0, b0)
    } else {
      log_ratio_current_ssr_matrixinv_2 <- log_ratio_approx(current_xi_2, eta_2, X, X_transpose, y, a0, b0, active_set_2)
      log_ratio_proposed_ssr_matrixinv_2 <- log_ratio_approx(proposed_xi_2, eta_2, X, X_transpose, y, a0, b0, active_set_2)
    }
    
    log_u <- log(runif(1))
    
    log_accept_prob_1 <- (log_ratio_proposed_ssr_matrixinv_1$log_likelihood - log_ratio_current_ssr_matrixinv_1$log_likelihood) + (log(proposed_xi_1)-log(current_xi_1))
    if (log_u<log_accept_prob_1){
      current_xi_1 <- proposed_xi_1
      log_ratio_current_ssr_matrixinv_1 <- log_ratio_proposed_ssr_matrixinv_1
    }
    
    log_accept_prob_2 <- (log_ratio_proposed_ssr_matrixinv_2$log_likelihood - log_ratio_current_ssr_matrixinv_2$log_likelihood) + (log(proposed_xi_2)-log(current_xi_2))
    if (log_u<log_accept_prob_2){
      current_xi_2 <- proposed_xi_2
      log_ratio_current_ssr_matrixinv_2 <- log_ratio_proposed_ssr_matrixinv_2
    }
  }
  return(list('xi_values'=c(current_xi_1, current_xi_2), 'log_ratio_ssr_matrix_inv_1'=log_ratio_current_ssr_matrixinv_1,  'log_ratio_ssr_matrix_inv_2'=log_ratio_current_ssr_matrixinv_2, 'active_set_1'= active_set_1, 'active_set_2'= active_set_2))
}

## Couplings for sigma^2 update given eta ##
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}
rigamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}
rigamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  x <- rigamma(1, alpha1, beta1)
  if (digamma(x, alpha1, beta1) + log(runif(1)) < digamma(x, alpha2, beta2)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rigamma(1, alpha2, beta2)
      reject <- (digamma(y, alpha2, beta2) + log(runif(1)) < digamma(y, alpha1, beta1))
    }
    return(c(x,y))
  }
}
### sigma^2 maximal coupling update step given eta and xi
sigma2_update_maximal_coupling <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0, sigma2_fixed_value=NULL)
{
  if(!is.null(sigma2_fixed_value)){
    sample <- c(sigma2_fixed_value, sigma2_fixed_value)
  } else {
    sample <- rigamma_coupled(((n+a0)/2), ((n+a0)/2), (ssr_1/2), (ssr_2/2))
  }
  return(sample)
}
### sigma^2 maximal coupling update step given eta and xi
sigma2_update_crn <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0,
                              sigma2_fixed_value=NULL)
{
  if(!is.null(sigma2_fixed_value)){
    sample <- c(1/sigma2_fixed_value, 1/sigma2_fixed_value)
  } else {
    # Common random number gamma draw
    crn_gamma <-  rgamma(1, shape = (n+a0)/2, rate = 1)
    sample <- crn_gamma / c((ssr_1/2), (ssr_2/2))
  }
  return(1/sample)
}
# Coupled update of sigma2 given eta
crn_max_sigma2_coupling <- function(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0,
                                    epsilon_sigma2=0, sigma2_fixed_value=NULL){
  relative_error_delta <- abs(ssr_1-ssr_2)
  if((0 < relative_error_delta) & (relative_error_delta < epsilon_sigma2)){
    output <- sigma2_update_maximal_coupling(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0, sigma2_fixed_value)
  } else {
    output <- sigma2_update_crn(xi_1, eta_1, xi_2, eta_2, n, ssr_1, ssr_2, a0, b0, sigma2_fixed_value)
  }
  return(output)
}

## Common random numbers coupling of beta given xi, eta, sigma2 ##
crn_joint_beta_update <-
  function(xi_1, sigma2_1, eta_1, xi_2, sigma2_2, eta_2,
           X, X_transpose, y, M_matrix_inverse_1, M_matrix_inverse_2,
           active_set_1, active_set_2)
  {
    n <- nrow(X)
    p <- nrow(X_transpose)
    # Using same common random numbers for draws on two chains
    random_u <- rnorm(p, 0, 1)
    random_delta <- c(rnorm(n,0,1))
    u_1 = (sqrt(xi_1*eta_1)^-1)*random_u
    v_1 = X%*%u_1 + random_delta
    v_star_1 <- M_matrix_inverse_1%*%(y/sqrt(sigma2_1) - v_1)
    if(sum(active_set_1)>0){
      U_1 = (xi_1^-1)*((eta_1[active_set_1]^(-1))*(X_transpose[active_set_1, ,drop=F]))
      u_1[active_set_1] <- u_1[active_set_1] + U_1%*%v_star_1
    }
    beta_parameter_1 <- sqrt(sigma2_1)*(u_1)
    u_2 = (sqrt(xi_2*eta_2)^-1)*random_u
    v_2 = X%*%u_2 + random_delta
    v_star_2 <- M_matrix_inverse_2%*%(y/sqrt(sigma2_2) - v_2)
    if(sum(active_set_2)>0){
      U_2 = (xi_2^-1)*((eta_2[active_set_2]^(-1))*(X_transpose[active_set_2, ,drop=F]))
      u_2[active_set_2] <- u_2[active_set_2] + U_2%*%v_star_2
    }
    beta_parameter_2 <- sqrt(sigma2_2)*(u_2)
    return(cbind(beta_parameter_1, beta_parameter_2))
  }

## Full coupled blocked Gibbs samplers ##
coupled_half_t_kernel <-
  function(X, X_transpose, y, a0=1, b0=1, std_MH=0.8,
           xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
           beta_1_current, beta_2_current, eta_1_current, eta_2_current,
           approximate_algo_delta_1=0, approximate_algo_delta_2=0, epsilon_eta = 0.5,
           epsilon_xi = Inf, epsilon_sigma2=Inf, nrepeats_eta=1,
           verbose = FALSE, xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df, metric)
  {
    n <- dim(X)[1]
    p <- dim(X)[2]
    
    if (verbose) ptm <- proc.time()
    # When epsilon_eta >=1, relative_error_delta <= epsilon_eta always.
    if(epsilon_eta >= 1){
      relative_error_delta <- 1
    } else {
      relative_error_delta <-
        half_t_max_couple_prob(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                               xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                               t_dist_df, iterations=1)
      # Checking overflow
      # if(typeof(relative_error_delta)=='list') return(relative_error_delta)
    }
    if (verbose) print(relative_error_delta)
    
    if (relative_error_delta <= epsilon_eta){ # Using max coupling of 1-step slice sampling when close
      eta_sample <-
        eta_update_half_t_max_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                     xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                     t_dist_df)
    } else {
      if (is.infinite(nrepeats_eta)){ stop("Number of slice sampling must be finite") }
      else {
        for (i in 1:nrepeats_eta) {
          eta_sample <-
            eta_update_half_t_crn_couple(xi_1_current, beta_1_current, eta_1_current, sigma2_1_current,
                                         xi_2_current, beta_2_current, eta_2_current, sigma2_2_current,
                                         t_dist_df)
          eta_1_current <- eta_sample[,1]
          eta_2_current <- eta_sample[,2]
        }
      }
    }
    eta_1_new <- eta_sample[,1]
    eta_2_new <- eta_sample[,2]
    if (verbose) print(proc.time()[3]-ptm[3])
    
    if (verbose) print(c('Eta components coupled:', sum(eta_1_new==eta_2_new)))
    
    if (verbose) ptm <- proc.time()
    xi_sample <-
      crn_max_xi_coupling(xi_1_current, eta_1_new, xi_2_current, eta_2_new,
                          X, X_transpose, y, a0, b0, std_MH,
                          approximate_algo_delta_1, approximate_algo_delta_2,
                          epsilon_xi, fixed=xi_fixed)
    xi_1_new <- xi_sample$xi_values[1]
    xi_2_new <- xi_sample$xi_values[2]
    if (verbose) print(proc.time()[3]-ptm[3])
    
    if (verbose) ptm <- proc.time()
    sigma2_fixed_value <- NULL
    if(sigma2_fixed==TRUE){sigma2_fixed_value <- sigma2_1_current}
    sigma2_sample <-
      crn_max_sigma2_coupling(xi_1_new,eta_1_new,xi_2_new,eta_2_new,n,
                              (xi_sample$log_ratio_ssr_matrix_inv_1)$ssr,
                              (xi_sample$log_ratio_ssr_matrix_inv_2)$ssr,a0,b0,
                              epsilon_sigma2,sigma2_fixed_value)
    sigma2_1_new <- sigma2_sample[1]
    sigma2_2_new <- sigma2_sample[2]
    if (verbose) print(proc.time()[3]-ptm[3])
    
    if (verbose) ptm <- proc.time()
    M_inverse_1 <- (xi_sample$log_ratio_ssr_matrix_inv_1)$M_matrix_inverse
    M_inverse_2 <- (xi_sample$log_ratio_ssr_matrix_inv_2)$M_matrix_inverse
    active_set_1 <- xi_sample$active_set_1
    active_set_2 <- xi_sample$active_set_2
    beta_samples <- crn_joint_beta_update(xi_1_new, sigma2_1_new, eta_1_new,
                                          xi_2_new, sigma2_2_new, eta_2_new,
                                          X, X_transpose, y, M_inverse_1, M_inverse_2,
                                          active_set_1, active_set_2)
    beta_1_new <- beta_samples[,1]
    beta_2_new <- beta_samples[,2]
    if (verbose) print(proc.time()[3]-ptm[3])
    
    chain1 <- list('beta'=beta_1_new, 'eta'=eta_1_new, 'sigma2'=sigma2_1_new, 'xi'=xi_1_new)
    chain2 <- list('beta'=beta_2_new, 'eta'=eta_2_new, 'sigma2'=sigma2_2_new, 'xi'=xi_2_new)
    
    metric_value <- metric(chain1, chain2)
    
    output <- list('beta_1_samples'=beta_1_new, 'beta_2_samples'=beta_2_new,
                   'eta_1_samples'=eta_1_new, 'eta_2_samples'=eta_2_new,
                   'sigma2_1_samples'=sigma2_1_new, 'sigma2_2_samples'=sigma2_2_new,
                   'xi_1_samples'=xi_1_new, 'xi_2_samples'=xi_2_new,
                   'metric'=metric_value)
    return(output)
  }


#' coupled_half_t_mcmc
coupled_half_t_mcmc <-
  function(mc_chain_size, X, X_transpose, y, a0=1, b0=1, std_MH=0.8, rinit=NULL,
           approximate_algo_delta_1=0, approximate_algo_delta_2=0, 
           epsilon_eta=0.5, epsilon_xi=Inf, epsilon_sigma2=Inf,
           nrepeats_eta=1, verbose = FALSE, totalduration = Inf, 
           xi_fixed=FALSE, sigma2_fixed=FALSE, t_dist_df, metric, store_states=FALSE){
    starttime <- Sys.time() # record starting time
    n <- dim(X)[1]
    p <- dim(X)[2]
    #
    if(is.null(rinit)){
      # Initializing from the prior
      rinit <- function(){
        xi <- (1/rt(1, df=1))^2
        sigma2 <- 1/rgamma(1, shape = a0/2, rate = b0/2)
        eta <- (1/rt(p, df=t_dist_df))^2
        beta <- rnorm(p)*sqrt(sigma2/(xi*eta))
        return(list(xi = xi, sigma2 = sigma2, beta = beta, eta = eta))
      }
    }
    #
    ## Initializing chains
    if(store_states){
      xi_samples1 <-     matrix(nrow = mc_chain_size, ncol = 1)
      sigma2_samples1 <- matrix(nrow = mc_chain_size, ncol = 1)
      beta_samples1 <-   matrix(nrow = mc_chain_size, ncol = p)
      eta_samples1 <-    matrix(nrow = mc_chain_size, ncol = p)
      xi_samples2 <-     matrix(nrow = mc_chain_size, ncol = 1)
      sigma2_samples2 <- matrix(nrow = mc_chain_size, ncol = 1)
      beta_samples2 <-   matrix(nrow = mc_chain_size, ncol = p)
      eta_samples2 <-    matrix(nrow = mc_chain_size, ncol = p)
    }
    metric_samples <- matrix(nrow = mc_chain_size, ncol = 1)
    #
    nrowsamples1 <- mc_chain_size
    # drawing initial states
    chain1 <- rinit()
    chain2 <- rinit()
    
    xi_1_current <-     chain1$xi
    sigma2_1_current <- chain1$sigma2
    beta_1_current <-   chain1$beta
    eta_1_current <-    chain1$eta
    xi_2_current <-     chain2$xi
    sigma2_2_current <- chain2$sigma2
    beta_2_current <-   chain2$beta
    eta_2_current <- chain2$eta
    
    if(store_states){
      xi_samples1[1,]  <-    xi_1_current
      sigma2_samples1[1,] <- sigma2_1_current
      beta_samples1[1,]  <-  beta_1_current
      eta_samples1[1,] <-    eta_1_current
      xi_samples2[1,]  <-    xi_2_current
      sigma2_samples2[1,] <- sigma2_2_current
      beta_samples2[1,]  <-  beta_2_current
      eta_samples2[1,] <- eta_2_current
    }
    
    chain1 <- list('beta'=beta_1_current, 'eta'=eta_1_current, 'sigma2'=sigma2_1_current, 'xi'=xi_1_current)
    chain2 <- list('beta'=beta_2_current, 'eta'=eta_2_current, 'sigma2'=sigma2_2_current, 'xi'=xi_2_current)
    metric_samples[1,] <- metric_current <- metric(chain1, chain2)
    
    # Setting up coupled chain
    iter <- 1
    meet <- FALSE
    finished <- FALSE
    meetingtime <- Inf
    while (!finished && iter < mc_chain_size){
      currentime <- Sys.time()
      elapsedtime <-
        as.numeric(lubridate::as.duration(lubridate::ymd_hms(currentime) -
                                            lubridate::ymd_hms(starttime)), "seconds")
      if (elapsedtime > totalduration){
        # time is up, interrupt function
        return(list(finished = FALSE, message = "interrupted because time is up"))
      }

      # res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, ...)
      output <-
        coupled_half_t_kernel(X, X_transpose, y, a0, b0, std_MH,
                              xi_1_current, xi_2_current, sigma2_1_current, sigma2_2_current,
                              beta_1_current, beta_2_current, eta_1_current, eta_2_current,
                              approximate_algo_delta_1 = approximate_algo_delta_1,
                              approximate_algo_delta_2 = approximate_algo_delta_2,
                              epsilon_eta = epsilon_eta, epsilon_xi = epsilon_xi, epsilon_sigma2 = epsilon_sigma2,
                              nrepeats_eta = nrepeats_eta, verbose = verbose,
                              xi_fixed=xi_fixed, sigma2_fixed=sigma2_fixed, t_dist_df, metric)
      
      xi_1_current <- output$xi_1_samples
      sigma2_1_current <- output$sigma2_1_samples
      beta_1_current <- output$beta_1_samples
      eta_1_current <- output$eta_1_samples
      xi_2_current <- output$xi_2_samples
      sigma2_2_current <- output$sigma2_2_samples
      beta_2_current <- output$beta_2_samples
      eta_2_current <- output$eta_2_samples
      
      if(store_states){
        xi_samples1[iter+1,] <- xi_1_current
        sigma2_samples1[iter+1,] <- sigma2_1_current
        beta_samples1[iter+1,] <- beta_1_current
        eta_samples1[iter+1,] <- eta_1_current
        xi_samples2[iter+1,] <- xi_2_current
        sigma2_samples2[iter+1,] <- sigma2_2_current
        beta_samples2[iter+1,] <- beta_2_current
        eta_samples2[iter+1,] <- eta_2_current
      }
      metric_samples[iter+1,] <- metric_current <- output$metric
      
      iter <- iter + 1
      if (iter >= mc_chain_size){finished <- TRUE}
      
      if(verbose) print(iter)
    }
    
    if(store_states){
      xi_samples1 <- xi_samples1[1:iter,,drop=F]
      xi_samples2 <- xi_samples2[1:iter,,drop=F]
      sigma2_samples1 <- sigma2_samples1[1:iter,,drop=F]
      sigma2_samples2 <- sigma2_samples2[1:iter,,drop=F]
      beta_samples1 <-  beta_samples1[1:iter,,drop=F]
      beta_samples2 <-  beta_samples2[1:iter,,drop=F]
      eta_samples1 <-   eta_samples1[1:iter,,drop=F]
      eta_samples2 <-   eta_samples2[1:iter,,drop=F]
    }
    metric_samples <- metric_samples[1:iter,,drop=F]
    
    if(store_states){
      final_output <- 
        list('beta_samples1'=beta_samples1, 'beta_samples2'=beta_samples2,
             'eta_samples1'=eta_samples1, 'eta_samples2'=eta_samples2,
             'sigma2_samples1'=sigma2_samples1, 'sigma2_samples2'=sigma2_samples2,
             'xi_samples1'=xi_samples1, 'xi_samples2'=xi_samples2,
             'metric_samples'=metric_samples, meetingtime = meetingtime, finished = TRUE)
    } else {
      final_output <- list('metric_samples'=metric_samples, meetingtime = meetingtime, finished = TRUE)
    }
    
    return(final_output)
  }



###############################################################################
#### Functions for the eta updates for blocked Gibbs sampler with half-t priors
###############################################################################
#### Single chain MCMC functions
## eta updates given beta, xi, sigma2
# Incomplete Gamma Function
low_inc_gamma <- function(rate, upper_truncation, log = TRUE){
  return(pgamma(upper_truncation,rate, log.p = log))
}
# Inverse of Lower incomplete Gamma Function
low_inc_gamma_inv <- function(rate, y, log = TRUE) {
  # If log=TRUE, the parameter y is taken to be on the natural logarithmic scale.
  if (log){
    if (any(y==0)){stop("Warnings: possible underflow")}
    if (any(y>0)){stop("Log second argument should be negative")}
  } else{
    if (any(y==1)){stop("Warnings: possible underflow")}
    if (any(y>1)){stop("Second argument should less than 1")}
  }
  return(qgamma(y, shape = rate, log.p = log))
}
# Confludent Hypergeomteric Function of the second kind function
conf_hyp <- function(a, b, z){
  integrand_ <- function(x) x^(a-1) * (1+x)^(b-a-1) * exp(-z*x)
  integrate(f = integrand_, lower = 0, upper = +Inf)$val / gamma(a)
}
# Lower Incomplete Confludent Hypergeometric Function of the second kind function
low_inc_conf_hyp <- function(a, b, z, t, regularize=TRUE){
  integrand_ <- function(x) x^(a-1) * (1+x)^(b-a-1) * exp(-z*x)
  if (regularize){
    integrate(f = integrand_, lower = 0, upper = t)$val /
      integrate(f = integrand_, lower = 0, upper = +Inf)$val
  } else{
    integrate(f = integrand_, lower = 0, upper = t)$val / gamma(a)
  }
}
# Perfect sampling from univariate p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]
r_trunc_poly_exp_crn <- function(poly_exponent, rate, trunc_upper, unif){
  u <- unif*(low_inc_gamma(poly_exponent, rate*trunc_upper, log = FALSE))
  return(low_inc_gamma_inv(poly_exponent, u, log=FALSE)/rate)
}
r_trunc_poly_exp <- function(poly_exponent, rate, trunc_upper){
  unif <- runif(length(rate))
  return(r_trunc_poly_exp_crn(poly_exponent, rate, trunc_upper, unif))
}
# Log Pdf of univariate p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]
d_trunc_poly_exp <- function(poly_exponent, rate, trunc_upper, eta){
  if(eta>trunc_upper) return(-Inf)
  log_d <- (poly_exponent-1)*log(eta)-rate*eta
  log_norm_constant <- -poly_exponent*log(rate)+lgamma(poly_exponent)+
    low_inc_gamma(poly_exponent, rate*trunc_upper, log = TRUE)
  return(log_d-log_norm_constant)
}
# Slice sampling from univariate p(x) \prop x^((v-1)/2)/(1+vx)^((v+1)/2)*exp(-mx) on [L,Inf]
eta_update_half_t <- function(xi, sigma2, beta, eta, t_dist_df, nrepeats=1)
{
  rate <- (beta^2)*(xi)/(2*sigma2)
  p <- length(eta)
  if (is.infinite(nrepeats)){
    stop('Perfect sampling not implemented for general t distribution shrinkage priors')
  } else {
    for (irepeat in 1:nrepeats){
      u <- runif(p)/(1 + t_dist_df*eta)^((1+t_dist_df)/2)
      eta <- r_trunc_poly_exp((1+t_dist_df)/2, rate,
                              (u^(-2/(1+t_dist_df))-1)/t_dist_df)
    }
  }
  return(eta)
}


###############################################################################
## CRN Coupling of Eta Update
eta_update_half_t_crn_couple <- function(xi_1, Beta_1, eta_1, sigma2_1,
                                         xi_2, Beta_2, eta_2, sigma2_2,
                                         t_dist_df){
  p <- length(eta_1)
  rate_1 <- (Beta_1^2)*(xi_1)/(2*sigma2_1)
  rate_2 <- (Beta_2^2)*(xi_2)/(2*sigma2_2)
  unif_crn_1 <- runif(p)
  u_1 <- unif_crn_1/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2)
  u_2 <- unif_crn_1/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2)
  unif_crn_2 <- runif(p)
  eta_1 <- r_trunc_poly_exp_crn((1+t_dist_df)/2, rate_1, (u_1^(-2/(1+t_dist_df))-1)/t_dist_df, unif_crn_2)
  eta_2 <- r_trunc_poly_exp_crn((1+t_dist_df)/2, rate_2, (u_2^(-2/(1+t_dist_df))-1)/t_dist_df, unif_crn_2)
  return(cbind(eta_1, eta_2))
}

## Maximal-coupling based Eta update
# Max couplings of univariates of form p(x) \prop x^(poly_exponent-1)*exp(-rate*x) on [0,trunc_upper]
trunc_poly_exp_max_couple <- function(poly_exponent, rate_1, rate_2,
                                      trunc_upper_1, trunc_upper_2){
  rp <- function() r_trunc_poly_exp(poly_exponent, rate_1, trunc_upper_1)
  rq <- function() r_trunc_poly_exp(poly_exponent, rate_2, trunc_upper_2)
  dp <- function(x){return(d_trunc_poly_exp(poly_exponent, rate_1, trunc_upper_1, x))}
  dq <- function(x){return(d_trunc_poly_exp(poly_exponent, rate_2, trunc_upper_2, x))}
  f <- get_max_coupling(rp, dp, rq, dq)
  return(f())
}
# Maximal-coupling based Eta update
eta_update_half_t_max_couple <- function(xi_1, Beta_1, eta_1, sigma2_1,
                                         xi_2, Beta_2, eta_2, sigma2_2,
                                         t_dist_df){
  p <- length(eta_1)
  rate_1 <- (Beta_1^2)*(xi_1)/(2*sigma2_1)
  rate_2 <- (Beta_2^2)*(xi_2)/(2*sigma2_2)
  etas_ <- matrix(0, nrow = p, ncol = 2)
  crn_unif <- runif(p)
  coupled_u <- cbind(crn_unif/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2),
                     crn_unif/(1 + t_dist_df*eta_2)^((1+t_dist_df)/2))
  truncs_ <- ((coupled_u)^(-2/(1+t_dist_df))-1)/t_dist_df
  for (j in 1:p){
    etas_[j,] <-
      trunc_poly_exp_max_couple((1+t_dist_df)/2, rate_1[j], rate_2[j],
                                truncs_[j,1], truncs_[j,2])
  }
  eta_1 <- etas_[,1]
  eta_2 <- etas_[,2]
  return(cbind(eta_1, eta_2))
}

## Distance metric
# Total variation between two univariates p(x) \prop x^(r-1)*exp(-mx) on [0,T]
# with r poly_exponent; m rate_1, rate_2; T trunc_upper_1, trunc_upper_2.
# Note: function is vectorized.
trunc_poly_exp_tv <- function(poly_exponent, rate_1, rate_2,
                              trunc_upper_1, trunc_upper_2){
  # Probabilities of coupling each component
  coupling_probs <- rep(NA, length(rate_1))
  
  rate_index1 <- rate_1==rate_2
  rate_index2 <- rate_1<rate_2
  rate_index3 <- rate_1>rate_2
  
  K <- (-poly_exponent*log(rate_1)+low_inc_gamma(poly_exponent, rate_1*trunc_upper_1, log = TRUE))-
    (-poly_exponent*log(rate_2)+low_inc_gamma(poly_exponent, rate_2*trunc_upper_2, log = TRUE))
  K[rate_2!=rate_1] <- (K/(rate_2-rate_1))[rate_2!=rate_1]
  trunc_upper_min <- exp(((log(trunc_upper_1) + log(trunc_upper_2)) -
                            abs(log(trunc_upper_1)-log(trunc_upper_2)))/2)
  trunc_upper_max <- exp(((log(trunc_upper_1) + log(trunc_upper_2)) +
                            abs(log(trunc_upper_1)-log(trunc_upper_2)))/2)
  
  K_index1 <- K<=0
  K_index2 <- K>=trunc_upper_min
  K_index3 <- (0<K)&(K<trunc_upper_min)
  
  # Calculating coupling probabilities case by case for K_index and rate_index
  coupling_probs[rate_index1] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[rate_2==rate_1], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_max)[rate_2==rate_1], log = TRUE))
  coupling_probs[K_index1 & rate_index2] <-
    exp(low_inc_gamma(poly_exponent,(rate_2*trunc_upper_min)[K_index1 & rate_index2], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index1 & rate_index2], log = TRUE))
  coupling_probs[K_index2 & rate_index2] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[K_index2 & rate_index2], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index2 & rate_index2], log = TRUE))
  coupling_probs[K_index3 & rate_index2] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*K)[K_index3 & rate_index2], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index3 & rate_index2], log = TRUE)) +
    exp(log(low_inc_gamma(poly_exponent,(rate_2*trunc_upper_min)[K_index3 & rate_index2], log = FALSE)-
              low_inc_gamma(poly_exponent,(rate_2*K)[K_index3 & rate_index2], log = FALSE))-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index3 & rate_index2], log = TRUE))
  coupling_probs[K_index1 & rate_index3] <-
    exp(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[K_index1 & rate_index3], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index1 & rate_index3], log = TRUE))
  coupling_probs[K_index2 & rate_index3] <-
    exp(low_inc_gamma(poly_exponent,(rate_2*trunc_upper_min)[K_index2 & rate_index3], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index2 & rate_index3], log = TRUE))
  coupling_probs[K_index3 & rate_index3] <-
    exp(low_inc_gamma(poly_exponent,(rate_2*K)[K_index3 & rate_index3], log = TRUE)-
          low_inc_gamma(poly_exponent,(rate_2*trunc_upper_2)[K_index3 & rate_index3], log = TRUE)) +
    exp(log(low_inc_gamma(poly_exponent,(rate_1*trunc_upper_min)[K_index3 & rate_index3], log = FALSE)-
              low_inc_gamma(poly_exponent,(rate_1*K)[K_index3 & rate_index3], log = FALSE))-
          low_inc_gamma(poly_exponent,(rate_1*trunc_upper_1)[K_index3 & rate_index3], log = TRUE))
  
  # Correcting rounding errors
  coupling_probs[coupling_probs<0] <- 0
  coupling_probs[coupling_probs>1] <- 1
  
  if(is.null(nrow(coupling_probs))){return(1-exp(sum(log(coupling_probs))))}
  return(1-exp(rowSums(log(coupling_probs))))
}
# TV UB from componentwise coupling
half_t_max_couple_prob <- function(xi_1, Beta_1, eta_1, sigma2_1,
                                   xi_2, Beta_2, eta_2, sigma2_2,
                                   t_dist_df, iterations=1){
  p <- length(eta_1)
  rate_1 <- (Beta_1^2)*(xi_1)/(2*sigma2_1)
  rate_2 <- (Beta_2^2)*(xi_2)/(2*sigma2_2)
  
  tv_ub <- rep(0, iterations)
  for (i in 1:iterations){
    crn_unif <- runif(p)
    coupled_u <- cbind(crn_unif/(1 + t_dist_df*eta_1)^((1+t_dist_df)/2),
                       crn_unif/(1 + t_dist_df*eta_2)^((1+t_dist_df)/2))
    truncs_ <- ((coupled_u)^(-2/(1+t_dist_df))-1)/t_dist_df
    tv_ub[i] <-
      trunc_poly_exp_tv((1+t_dist_df)/2, rate_1, rate_2, truncs_[,1], truncs_[,2])
    
    if(is.na(tv_ub[i])) # Checking underflow
    {
      print(c("Warning: possible overflow", i))
      return(list(xi_1=xi_1, Beta_1=Beta_1, eta_1=eta_1, sigma2_1=sigma2_1,
                  xi_2=xi_2, Beta_2=Beta_2, eta_2=eta_2, sigma2_2=sigma2_2,
                  p=p, eta_lower_bound=eta_lower_bound, t_dist_df=t_dist_df,
                  iterations=iterations, coupled_u=coupled_u))
    }
  }
  return(mean(tv_ub))
}

###############################################################################
#### Functions for the maximal couplings
###############################################################################
get_max_coupling <- function(rp, dp, rq, dq){
  function(){
    x <- rp()
    logu <- log(runif(1))
    if (dp(x) + logu < dq(x)){
      return(c(x,x))
    } else {
      reject <- TRUE
      y <- NA
      while (reject){
        y <- rq()
        logu <- log(runif(1))
        reject <- (dq(y) + logu < dp(y))
      }
      return(c(x,y))
    }
  }
}

### sample from maximal coupling of truncated exponentials (truncated on [0,trunc])
### with rate rate1, rate2, and truncations trunc1, trunc2
rtruncexp_maxcoupling <- function(rate1, trunc1, rate2, trunc2, eta_lower_bound){
  rp <- function() gen_truncated_exp(rate1, trunc1, eta_lower_bound)
  rq <- function() gen_truncated_exp(rate2, trunc2, eta_lower_bound)
  dp <- function(x){
    if (x > trunc1) return(-Inf) 
    else return(dexp((x-eta_lower_bound), rate1, log = TRUE) - pexp((trunc1-eta_lower_bound), rate1, log.p = TRUE))
  }
  dq <- function(x){
    if (x > trunc2) return(-Inf)
    else return(dexp((x-eta_lower_bound), rate2, log = TRUE) - pexp((trunc2-eta_lower_bound), rate2, log.p = TRUE))
  }
  f <- get_max_coupling(rp, dp, rq, dq)
  return(f())
}

### Sample from maximal coupling of two uniform distribution on [0, s] and [0,t]
### For s<t, sample U~Unif(0,s) and set X=sU. Set Y=sU w.p. s/t, Y=(t-s)U+s w.p. 1-s/t
runif_maxcoupling <- function(s, t){
  trun_min <- exp((log(s)+log(t) - abs(log(s)-log(t)))/2) # calculate min(s,t)
  trun_min[(s+t)==0] <- 0
  trun_max <- exp((log(s)+log(t) + abs(log(s)-log(t)))/2) # calculate max(s,t)
  p <- length(s)
  unif_draw <- runif(p)
  
  x <- unif_draw*trun_min
  y <- x
  coins <- runif(p)
  y[coins>trun_min/trun_max] <-
    (trun_max*unif_draw-x+trun_min)[coins>trun_min/trun_max]
  
  # Swapping back for all components
  swap_bool <- s>t
  x[swap_bool] <- x[swap_bool] + y[swap_bool]
  y[swap_bool] <- x[swap_bool] - y[swap_bool]
  x[swap_bool] <- x[swap_bool] - y[swap_bool]
  return(cbind(x,y))
}
