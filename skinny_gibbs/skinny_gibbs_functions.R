# Functions for Skinny Gibbs implementation

require(Rcpp)
require(RcppEigen)
# Rcpp functions
Rcpp::cppFunction('
Eigen::MatrixXd cpp_mat_vec_prod(const Eigen::MatrixXd X, const Eigen::VectorXd y){
Eigen::VectorXd output;
output = X*y;
return output;
}', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd fcprd(const Eigen::MatrixXd X){
const int n = X.cols();
  return Eigen::MatrixXd(n, n).setZero().selfadjointView<Eigen::Lower>().rankUpdate(X.adjoint());
}', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd cpp_prod(const Eigen::MatrixXd X, const Eigen::MatrixXd Y){
  return Eigen::MatrixXd(X*Y);
}', depends = 'RcppEigen')
Rcpp::cppFunction('
Eigen::MatrixXd cpp_dotprod(const Eigen::VectorXd X, const Eigen::VectorXd Y){
  return Eigen::VectorXd(X.transpose()*Y);
}', depends = 'RcppEigen')



# Exact chain functions
# beta update
update_beta_exact <- function(w, z, aug_y, X, tau0, tau1, u=NULL, delta=NULL){
  d <- as.vector(z/(tau1^2)+(1-z)/(tau0^2))
  # wX <- X*(w^0.5)
  # wX_transpose <- t(wX)

  w <- as.vector(w)
  wX_transpose <- t(X)*(w^0.5)
  wX <- t(wX_transpose)
  
  w_aug_y <- (aug_y)*(w^0.5)
  
  p <- length(z)
  n <- length(aug_y)
  
  if(is.null(u)){u <- rnorm(p, 0, 1)}
  u = u/sqrt(d)
  if(is.null(delta)){delta <- c(rnorm(n,0,1))}
  v = wX%*%u + delta
  
  wXd_transpose <- wX_transpose/d^0.5
  M_matrix <- fcprd(wXd_transpose)+diag(n)
  Mchol <- chol(M_matrix)
  M_matrix_inverse <- chol2inv(Mchol)
  
  v_star <- M_matrix_inverse%*%(w_aug_y - v)
  beta = u + (d^(-1))*cpp_mat_vec_prod(wX_transpose,v_star)
  return(beta)
}
# aug_y update
update_aug_y_exact <- function(beta, w, y, X, u_crn=NULL){
  n <- length(y)
  if(is.null(u_crn)){u_crn <- runif(n)}
  means <- cpp_mat_vec_prod(X, beta)
  sds <- 1/sqrt(w)
  ubs <- rep(Inf,n)
  lbs <- rep(-Inf,n)
  lbs[y==1] <- 0
  ubs[y!=1] <- 0
  aug_y <- TruncatedNormal::qtnorm(u_crn, mu=means, sd=sds, lb=lbs, ub=ubs)
  return(aug_y)
}
# z update
update_z_exact <- function(beta, w, tau0, tau1, q, u_crn=NULL){
  p <- length(beta)
  if(is.null(u_crn)){u_crn <- runif(p)}
  prob1 <- q*dnorm(beta, sd=tau1)
  prob2 <- (1-q)*dnorm(beta, sd=tau0)
  probs <- prob1/(prob1+prob2)
  z <- ifelse(u_crn<probs,1,0)
  return(z)
}
# w update
update_w_exact <- function(beta, X, y, nu, u_crn=NULL){
  n <- length(y)
  if(is.null(u_crn)){u_crn <- runif(n)}
  wsquared <- (pi^2)*(nu-2)/(3*nu)
  residuals <- y-cpp_mat_vec_prod(X, beta)
  gamma_samples <- qgamma(u_crn, shape=(nu+1)/2, rate=1)
  return(gamma_samples/((wsquared*nu+residuals^2)/2))
}



# Exact chain kernel with w fixed
full_gibbs_kernel <- function(beta, aug_y, z, w, X, y, tau0, tau1, q, nu, 
                              w_fixed, random_samples){
  beta_new <- 
    update_beta_exact(w, z, aug_y, X, tau0, tau1, u=random_samples$beta_u,
                      delta=random_samples$beta_delta)
  aug_y_new <- update_aug_y_exact(beta_new, w, y, X, u_crn=random_samples$aug_y_u)
  z_new <- update_z_exact(beta_new, w, tau0, tau1, q, u_crn=random_samples$z_u)
  
  if(w_fixed){
    w_new <- w
  } else {
    w_new <- update_w_exact(beta_new, X, y, nu, u_crn=random_samples$w_u)
  }
  return(list('beta'=beta_new, 'aug_y'=aug_y_new, 'z'=z_new, 'w'=w_new))
}

# Exact chain
full_gibbs <- function(X, y, tau0, tau1, q, nu=7.3, chain_length=1e3, rinit=NULL,
                       w_fixed=TRUE){
  p <- dim(X)[2]
  n <- length(y)
  if(is.null(rinit)){
    # Sampling from the prior
    z <- rbinom(p,1,q)
    beta <- rnorm(p)
    beta[z==0] <- beta[z==0]*tau0
    beta[z==1] <- beta[z==1]*tau1
    
    if(w_fixed){ # Fixing w
      s <- pi/sqrt(3) # sd of logistic distribution
      w <- rep(1,n)/(s^2)
    } else {
      w <- update_w_exact(beta, X, y, nu)
    }
    aug_y <- X%*%beta + w*rnorm(n, mean = 0, sd = 1)
  }
  beta_samples <- matrix(NA, nrow = chain_length, ncol = p)
  aug_y_samples <- matrix(NA, nrow = chain_length, ncol = n)
  z_samples <- matrix(NA, nrow = chain_length, ncol = p)
  w_samples <- matrix(NA, nrow = chain_length, ncol = n)
  for(t in 1:chain_length){
    random_samples <- 
      list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), 
           aug_y_u=runif(n), z_u=runif(p), w_u=runif(n))
    new_state <- 
      full_gibbs_kernel(beta, aug_y, z, w, X, y, tau0, tau1, q, nu, w_fixed, random_samples)
    beta <- new_state$beta
    aug_y <- new_state$aug_y
    z <- new_state$z
    w <- new_state$w
    
    beta_samples[t,] <- beta
    aug_y_samples[t,] <- aug_y
    z_samples[t,] <- z
    w_samples[t,] <- w
  }
  return(list('beta'=beta_samples, 'aug_y'=aug_y_samples, 'z'=z_samples, 'w'=w_samples))
}

# Skinny chain functions
# beta update
update_beta_skinny <- function(w, z, aug_y, X, tau0, tau1, u=NULL, delta=NULL,
                               X_scaled=TRUE){
  p <- length(z)
  n <- length(aug_y)
  if(is.null(u)){u <- rnorm(p, 0, 1)}
  if(is.null(delta)){delta <- c(rnorm(n,0,1))}
  
  active_set <- z==1
  
  if(sum(active_set)==p){
    beta_out <- update_beta_exact(w, z, aug_y, X, tau0, tau1, u, delta)
  } else if (sum(active_set)==0){
    if(X_scaled){
      beta_out <- u/sqrt((n-1) + tau0^(-2))
    } else{
      X_diags <- sapply(c(1:p), function(x) sum(X[,x]^2))
      beta_out <- u/sqrt(X_diags + tau0^(-2))
    }
  } else {
    z_active <- z[active_set]
    X_active <- X[,active_set,drop=FALSE]
    u_active <- u[active_set]
    beta_active <- update_beta_exact(w,z_active,aug_y,X_active, 
                                     tau0, tau1, u_active, delta)
    
    u_inactive <- u[!active_set]
    X_inactive <- X[,!active_set,drop=FALSE]
    if(X_scaled){
      beta_inactive <- u_inactive/sqrt((n-1) + tau0^(-2))
    } else{
      X_diags <- sapply(c(1:sum(!active_set)), function(x) sum(X_inactive[,x]^2))
      beta_inactive <- u_inactive/sqrt(X_diags + tau0^(-2))
    }
    
    beta_out <- rep(NA,p)
    beta_out[active_set] <- beta_active
    beta_out[!active_set] <- beta_inactive
  }
  return(beta_out)
}
# aug_y update
update_aug_y_skinny <- function(beta, w, z, y, X, u_crn=NULL){
  n <- length(y)
  if(is.null(u_crn)){u_crn <- runif(n)}
  
  active_set <- z==1
  
  means <- cpp_mat_vec_prod(X[,active_set, drop=FALSE], beta[active_set])
  sds <- 1/sqrt(w)
  ubs <- rep(Inf,n)
  lbs <- rep(-Inf,n)
  lbs[y==1] <- 0
  ubs[y!=1] <- 0
  aug_y <- TruncatedNormal::qtnorm(u_crn, mu=means, sd=sds, lb=lbs, ub=ubs)
  return(aug_y)
}
# z update
update_z_skinny <- function(beta, w, z, aug_y, X, tau0, tau1, q, 
                            X_scaled, w_fixed, u_crn=NULL){
  beta <- as.vector(beta)
  
  p <- length(beta)
  if(is.null(u_crn)){u_crn <- runif(p)}
  
  log_prob1 <- log(q)+dnorm(beta, sd=tau1, log = TRUE)
  log_prob0 <- log(1-q)+dnorm(beta, sd=tau0, log=TRUE)
  
  active <- z==1
  X_active <- X[,active,drop=FALSE]
  beta_active <- beta[active]
  residuals <- aug_y-cpp_mat_vec_prod(X_active, beta_active)
  beta_X_colunnwise_product <- t(t(X)*beta)
  for (j in 1:p){
    # Xjbetaj <- beta[j]*X[,j]
    Xjbetaj <- beta_X_colunnwise_product[,j]
    residuals_j <- residuals + Xjbetaj
    
    # if(X_scaled){
    #   modification <- # ONLY when w has one unique value and X is scaled
    #     sum(Xjbetaj*residuals_j)*w[1] + 0.5*beta[j]^2*(n-1)*(1-w[1])
    # } else{
    #   modification <-
    #     sum(Xjbetaj*w*residuals_j) + 0.5*beta[j]^2*sum(X[,j]*(1-w)*X[,j])
    # }

    # if(X_scaled){
    #   # ONLY when w has one unique value and X is scaled
    #   modification <-
    #     (cpp_dotprod(Xjbetaj,residuals_j)*w[1] + 0.5*beta[j]^2*(n-1)*(1-w[1]))[1,1]
    # } else{
    #   modification <-
    #     (cpp_dotprod(Xjbetaj,w*residuals_j) + 0.5*beta[j]^2*sum(X[,j]*(1-w)*X[,j]))[1,1]
    # }
    # 
    
    if(w_fixed){
      if(X_scaled){
        # When w has one unique value and X is scaled
        modification <-
          (cpp_dotprod(Xjbetaj,residuals_j)*w[1] + 0.5*beta[j]^2*(n-1)*(1-w[1]))[1,1]
      } else {
        modification <-
          (cpp_dotprod(Xjbetaj,residuals_j)*w[1] + 
             0.5*beta[j]^2*cpp_dotprod(X[,j],X[,j])*(1-w[1]))[1,1]
        # 0.5*beta[j]^2*sum(X[,j]*X[,j])*(1-w[1]))[1,1]
      }
    } else {
      modification <-
        (cpp_dotprod(Xjbetaj,w*residuals_j) + 
           0.5*beta[j]^2*cpp_dotprod(X[,j],(1-w)*X[,j]))[1,1]
      # 0.5*beta[j]^2*sum(X[,j]*(1-w)*X[,j]))[1,1]
    }
    
    modified_log_odds <- log_prob1[j] - log_prob0[j] + modification
    log_prob <- -log(1+exp(-modified_log_odds))
    if(log(u_crn[j])<log_prob){
      z[j] <- 1
      active[j] <- TRUE
      X_active <- X[,active,drop=FALSE]
      beta_active <- beta[active]
      residuals <- aug_y-cpp_mat_vec_prod(X_active, beta_active)
    } else{
      z[j] <- 0
    }
    # print(c(exp(log_prob),sum(active_minus_j),j))
  }
  return(z)
}
# w update
update_w_skinny <- function(beta, z, X, y, nu, u_crn=NULL){
  n <- length(y)
  if(is.null(u_crn)){u_crn <- runif(n)}
  
  active_set <- z==1
  
  wsquared <- (pi^2)*(nu-2)/(3*nu)
  residuals <- y-cpp_mat_vec_prod(X[,active_set, drop=FALSE], beta[active_set])
  gamma_samples <- qgamma(u_crn, shape=(nu+1)/2, rate=1)
  return(gamma_samples/((wsquared*nu+residuals^2)/2))
}

# Skinny chain kernel with w fixed
skinny_gibbs_kernel <- function(beta, aug_y, z, w, X, y, tau0, tau1, q, nu,
                                X_scaled, skinny_correction, w_fixed, random_samples){
  beta_new <- 
    update_beta_skinny(w, z, aug_y, X, tau0, tau1, X_scaled=X_scaled, 
                       u=random_samples$beta_u, delta=random_samples$beta_delta)
  aug_y_new <- update_aug_y_skinny(beta_new, w, z, y, X, u_crn=random_samples$aug_y_u)
  if(skinny_correction){ # Correction
    z_new <- update_z_skinny(beta_new, w, z, aug_y, X, tau0, tau1, q, 
                             X_scaled=X_scaled, w_fixed=w_fixed, u_crn=random_samples$z_u)
  } else { # No correction
    z_new <- update_z_exact(beta_new, w, tau0, tau1, q, u_crn=random_samples$z_u)
  }
  
  if(w_fixed){
    w_new <- w
  } else {
    w_new <- update_w_skinny(beta_new, z_new, X, y, nu, u_crn=random_samples$w_u)
  }
  
  return(list('beta'=beta_new, 'aug_y'=aug_y_new, 'z'=z_new, 'w'=w_new))
}

skinny_gibbs <- function(X, y, tau0, tau1, q, nu=7.3, chain_length=1e3, rinit=NULL,
                         X_scaled=TRUE, skinny_correction=TRUE, w_fixed=TRUE){
  p <- dim(X)[2]
  n <- length(y)
  if(is.null(rinit)){
    # Sampling from the prior
    z <- rbinom(p,1,q)
    beta <- rnorm(p)
    beta[z==0] <- beta[z==0]*tau0
    beta[z==1] <- beta[z==1]*tau1
    
    if(w_fixed){ # Fixing w
      s <- pi/sqrt(3) # sd of logistic distribution
      w <- rep(1,n)/(s^2)
    } else {
      w <- update_w_skinny(beta, z, X, y, nu)
    }
    
    aug_y <- X%*%beta + w*rnorm(n, mean = 0, sd = 1)
  }
  
  beta_samples <- matrix(NA, nrow = chain_length, ncol = p)
  aug_y_samples <- matrix(NA, nrow = chain_length, ncol = n)
  z_samples <- matrix(NA, nrow = chain_length, ncol = p)
  w_samples <- matrix(NA, nrow = chain_length, ncol = n)
  
  for(t in 1:chain_length){
    random_samples <- 
      list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), aug_y_u=runif(n), z_u=runif(p),
           w_u=runif(n))
    new_state <- 
      skinny_gibbs_kernel(beta, aug_y, z, w, X, y, tau0, tau1, q, nu,
                          X_scaled, skinny_correction, w_fixed, random_samples)
    beta <- new_state$beta
    aug_y <- new_state$aug_y
    z <- new_state$z
    w <- new_state$w
    
    beta_samples[t,] <- beta
    aug_y_samples[t,] <- aug_y
    z_samples[t,] <- z
    w_samples[t,] <- w
  }
  
  return(list('beta'=beta_samples, 'aug_y'=aug_y_samples, 'z'=z_samples, 
              'w'=w_samples))
}



# CRN coupled chain kernel
exact_skinny_coupled_kernel <- function(exact_beta, exact_aug_y, exact_z, exact_w, 
                                        skinny_beta, skinny_aug_y, skinny_z, skinny_w, 
                                        X, y, tau0, tau1, q, nu, X_scaled, 
                                        skinny_correction, w_fixed, random_samples){
  p <- length(exact_beta)
  n <- length(y)
  
  exact_sample_new <- 
    full_gibbs_kernel(exact_beta, exact_aug_y, exact_z, exact_w, X, y, 
                      tau0, tau1, q, nu, w_fixed, random_samples)
  skinny_sample_new <- 
    skinny_gibbs_kernel(skinny_beta, skinny_aug_y, skinny_z, skinny_w, X, y, 
                      tau0, tau1, q, nu, X_scaled, skinny_correction, w_fixed, random_samples)
  
  # u_beta <- rnorm(p)
  # delta_beta <- rnorm(n)
  # exact_beta_new <-
  #   update_beta_exact(exact_w, exact_z, exact_aug_y, X, tau0, tau1,
  #                     u=u_beta, delta=delta_beta)
  # skinny_beta_new <- 
  #   update_beta_skinny(skinny_w, skinny_z, skinny_aug_y, X, tau0, tau1, u=u_beta, 
  #                      delta=delta_beta, X_scaled=X_scaled)
  # u_aug_y <- runif(n)
  # exact_aug_y_new <- update_aug_y(exact_beta_new, exact_w, y, X, u_crn=u_aug_y)
  # skinny_aug_y_new <- update_aug_y_skinny(skinny_beta_new, skinny_w, skinny_z, y, X, u_crn=u_aug_y)
  # u_z <- runif(p)
  # exact_z_new <- update_z_exact(exact_beta_new, exact_w, tau0, tau1, q, u_crn=u_z)
  # if(skinny_correction){ # Correction
  #   skinny_z_new <- 
  #     update_z_skinny(skinny_beta_new, skinny_w, skinny_z, skinny_aug_y_new,
  #                     X, tau0, tau1, q, u_crn=u_z, X_scaled=X_scaled)
  # } else {  # No correction
  #   skinny_z_new <- 
  #     update_z_exact(skinny_beta_new, skinny_w, tau0, tau1, q, u_crn=u_z)
  # }
  # 
  # exact_w_new <- exact_w
  # skinny_w_new <- skinny_w
  
  return(list('exact_beta'=exact_sample_new$beta, 
              'exact_aug_y'=exact_sample_new$aug_y, 
              'exact_z'=exact_sample_new$z, 'exact_w'=exact_sample_new$w,
              'skinny_beta'=skinny_sample_new$beta, 
              'skinny_aug_y'=skinny_sample_new$aug_y,
              'skinny_z'=skinny_sample_new$z, 'skinny_w'=skinny_sample_new$w))
}

# CRN coupled chains
coupled_exact_chain <- 
  function(X, y, tau0, tau1, q, nu=7.3, chain_length=1e3, rinit=NULL, w_fixed=TRUE){
  p <- dim(X)[2]
  n <- length(y)
  if(is.null(rinit)){
    # Sampling from an initial distribution close to the prior
    # z1 <- rbinom(p,1,q)
    # Sampling from an initial distribution close with all z1 1s
    z1 <- rbinom(p,1,1)
    
    beta1 <- rnorm(p)
    beta1[z1==0] <- beta1[z1==0]*tau0
    beta1[z1==1] <- beta1[z1==1]*tau1
    if(w_fixed){ # Fixing w
      s <- pi/sqrt(3) # sd of logistic distribution
      w1 <- rep(1,n)/(s^2)
    } else {
      w1 <- update_w_exact(beta1, X, y, nu)
    }
    aug_y1 <- X%*%beta1 + w1*rnorm(n, mean = 0, sd = 1)
    
    # Sampling from an initial distribution close to the prior
    # z2 <- rbinom(p,1,q)
    # Sampling from an initial distribution close with all z1 0s
    z2 <- rbinom(p,1,0)
    
    beta2 <- rnorm(p)
    beta2[z2==0] <- beta2[z2==0]*tau0
    beta2[z2==1] <- beta2[z2==1]*tau1
    if(w_fixed){ # Fixing w
      s <- pi/sqrt(3) # sd of logistic distribution
      w2 <- rep(1,n)/(s^2)
    } else {
      w2 <- update_w_exact(beta2, X, y, nu)
    }
    aug_y2 <- X%*%beta2 + w2*rnorm(n, mean = 0, sd = 1)
  }
  
  beta_samples1 <- matrix(NA, nrow = chain_length, ncol = p)
  aug_y_samples1 <- matrix(NA, nrow = chain_length, ncol = n)
  z_samples1 <- matrix(NA, nrow = chain_length, ncol = p)
  w_sample1 <- matrix(NA, nrow = chain_length, ncol = n)
  
  beta_samples2 <- matrix(NA, nrow = chain_length, ncol = p)
  aug_y_samples2 <- matrix(NA, nrow = chain_length, ncol = n)
  z_samples2 <- matrix(NA, nrow = chain_length, ncol = p)
  w_sample2 <- matrix(NA, nrow = chain_length, ncol = n)
  
  for(t in 1:chain_length){
    random_samples <- 
      list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), aug_y_u=runif(n), 
           z_u=runif(p), w_u=runif(n))
    new_state1 <- full_gibbs_kernel(beta1, aug_y1, z1, w1, X, y, tau0, tau1, q, 
                                    nu, w_fixed, random_samples)
    new_state2 <- full_gibbs_kernel(beta2, aug_y2, z2, w2, X, y, tau0, tau1, q, 
                                    nu, w_fixed, random_samples)
    
    beta1 <- new_state1$beta
    aug_y1 <- new_state1$aug_y
    z1 <- new_state1$z
    w1 <- new_state1$w
    beta2 <- new_state2$beta
    aug_y2 <- new_state2$aug_y
    z2 <- new_state2$z
    w2 <- new_state2$w
    
    beta_samples1[t,] <- beta1
    aug_y_samples1[t,] <- aug_y1
    z_samples1[t,] <- z1
    w_sample1[t,] <- w1
    beta_samples2[t,] <- beta2
    aug_y_samples2[t,] <- aug_y2
    z_samples2[t,] <- z2
    w_sample2[t,] <- w2
  }
  return(list('beta1'=beta_samples1, 'aug_y1'=aug_y_samples1, 'z1'=z_samples1, 'w1'=w_sample1,
              'beta2'=beta_samples2, 'aug_y2'=aug_y_samples2, 'z2'=z_samples2, 'w2'=w_sample2))
}

# CRN coupled chains
coupled_skinny_chain <- 
  function(X, y, tau0, tau1, q, nu=7.3, chain_length=1e3, 
           rinit=NULL, X_scaled=TRUE, skinny_correction=TRUE, w_fixed=TRUE){
    p <- dim(X)[2]
    n <- length(y)
    if(is.null(rinit)){
      # Sampling from an initial distribution close to the prior
      # z1 <- rbinom(p,1,q)
      # Sampling from an initial distribution close with all z1 1s
      z1 <- rbinom(p,1,1)
      
      beta1 <- rnorm(p)
      beta1[z1==0] <- beta1[z1==0]*tau0
      beta1[z1==1] <- beta1[z1==1]*tau1
      if(w_fixed){ # Fixing w
        s <- pi/sqrt(3) # sd of logistic distribution
        w1 <- rep(1,n)/(s^2)
      } else {
        w1 <- update_w_skinny(beta1, z1, X, y, nu)
      }
      aug_y1 <- X%*%beta1 + w1*rnorm(n, mean = 0, sd = 1)
      
      # Sampling from an initial distribution close to the prior
      # z2 <- rbinom(p,1,q)
      # Sampling from an initial distribution close with all z1 1s
      z2 <- rbinom(p,1,0)
      
      beta2 <- rnorm(p)
      beta2[z2==0] <- beta2[z2==0]*tau0
      beta2[z2==1] <- beta2[z2==1]*tau1
      if(w_fixed){ # Fixing w
        s <- pi/sqrt(3) # sd of logistic distribution
        w2 <- rep(1,n)/(s^2)
      } else {
        w2 <- update_w_skinny(beta2, z2, X, y, nu)
      }
      aug_y2 <- X%*%beta2 + w2*rnorm(n, mean = 0, sd = 1)
    }
    
    beta_samples1 <- matrix(NA, nrow = chain_length, ncol = p)
    aug_y_samples1 <- matrix(NA, nrow = chain_length, ncol = n)
    z_samples1 <- matrix(NA, nrow = chain_length, ncol = p)
    w_sample1 <- matrix(NA, nrow = chain_length, ncol = n)
    
    beta_samples2 <- matrix(NA, nrow = chain_length, ncol = p)
    aug_y_samples2 <- matrix(NA, nrow = chain_length, ncol = n)
    z_samples2 <- matrix(NA, nrow = chain_length, ncol = p)
    w_sample2 <- matrix(NA, nrow = chain_length, ncol = n)
    
    for(t in 1:chain_length){
      random_samples <- 
        list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), aug_y_u=runif(n), 
             z_u=runif(p), w_u=runif(n))
      new_state1 <- 
        skinny_gibbs_kernel(beta1, aug_y1, z1, w1, X, y, tau0, tau1, q, nu,
                            X_scaled, skinny_correction, w_fixed, random_samples)
      new_state2 <- 
        skinny_gibbs_kernel(beta2, aug_y2, z2, w2, X, y, tau0, tau1, q, nu,
                            X_scaled, skinny_correction, w_fixed, random_samples)
      
      beta1 <- new_state1$beta
      aug_y1 <- new_state1$aug_y
      z1 <- new_state1$z
      w1 <- new_state1$w
      beta2 <- new_state2$beta
      aug_y2 <- new_state2$aug_y
      z2 <- new_state2$z
      w2 <- new_state2$w
      
      beta_samples1[t,] <- beta1
      aug_y_samples1[t,] <- aug_y1
      z_samples1[t,] <- z1
      w_sample1[t,] <- w1
      beta_samples2[t,] <- beta2
      aug_y_samples2[t,] <- aug_y2
      z_samples2[t,] <- z2
      w_sample2[t,] <- w2
    }
    return(list('beta1'=beta_samples1, 'aug_y1'=aug_y_samples1, 'z1'=z_samples1, 'w1'=w_sample1,
                'beta2'=beta_samples2, 'aug_y2'=aug_y_samples2, 'z2'=z_samples2, 'w2'=w_sample2))
  }

# CRN of exact and skinny
exact_skinny_chain <- function(X, y, tau0, tau1, q, nu=7.3, chain_length=1e3, rinit=NULL,
                               X_scaled=TRUE, skinny_correction=TRUE, w_fixed=TRUE){
  p <- dim(X)[2]
  n <- length(y)
  if(is.null(rinit)){
    # Sampling from an initial distribution close to the prior
    exact_z <- rbinom(p,1,q)
    exact_beta <- rnorm(p)
    exact_beta[exact_z==0] <- exact_beta[exact_z==0]*tau0
    exact_beta[exact_z==1] <- exact_beta[exact_z==1]*tau1
    if(w_fixed){ # Fixing w
      s <- pi/sqrt(3) # sd of logistic distribution
      exact_w <- rep(1,n)/(s^2)
    } else {
      exact_w <- update_w_exact(exact_beta, X, y, nu)
    }
    exact_aug_y <- X%*%exact_beta + exact_w*rnorm(n, mean = 0, sd = 1)
    
    skinny_z <- rbinom(p,1,q)
    skinny_beta <- rnorm(p)
    skinny_beta[skinny_z==0] <- skinny_beta[skinny_z==0]*tau0
    skinny_beta[skinny_z==1] <- skinny_beta[skinny_z==1]*tau1
    if(w_fixed){ # Fixing w
      s <- pi/sqrt(3) # sd of logistic distribution
      skinny_w <- rep(1,n)/(s^2)
    } else {
      skinny_w <- update_w_skinny(skinny_beta, skinny_z, X, y, nu)
    }
    skinny_aug_y <- X%*%exact_beta + skinny_w*rnorm(n, mean = 0, sd = 1)
  }
  
  exact_beta_samples <- matrix(NA, nrow = chain_length, ncol = p)
  exact_aug_y_samples <- matrix(NA, nrow = chain_length, ncol = n)
  exact_z_samples <- matrix(NA, nrow = chain_length, ncol = p)
  exact_w_samples <- matrix(NA, nrow = chain_length, ncol = n)
  
  skinny_beta_samples <- matrix(NA, nrow = chain_length, ncol = p)
  skinny_aug_y_samples <- matrix(NA, nrow = chain_length, ncol = n)
  skinny_z_samples <- matrix(NA, nrow = chain_length, ncol = p)
  skinny_w_samples <- matrix(NA, nrow = chain_length, ncol = n)
  
  for(t in 1:chain_length){
    random_samples <- 
      list(beta_u=rnorm(p, 0, 1), beta_delta=rnorm(n,0,1), aug_y_u=runif(n), z_u=runif(p))
    new_state <- 
      exact_skinny_coupled_kernel(exact_beta, exact_aug_y, exact_z, exact_w, 
                                  skinny_beta, skinny_aug_y, skinny_z, skinny_w, 
                                  X, y, tau0, tau1, q, nu, X_scaled, skinny_correction,
                                  w_fixed, random_samples)
    exact_beta <- new_state$exact_beta
    exact_aug_y <- new_state$exact_aug_y
    exact_z <- new_state$exact_z
    exact_w <- new_state$exact_w
    
    skinny_beta <- new_state$skinny_beta
    skinny_aug_y <- new_state$skinny_aug_y
    skinny_z <- new_state$skinny_z
    skinny_w <- new_state$skinny_w
    
    exact_beta_samples[t,] <- exact_beta
    exact_aug_y_samples[t,] <- exact_aug_y
    exact_z_samples[t,] <- exact_z
    exact_w_samples[t,] <- exact_w
    
    skinny_beta_samples[t,] <- skinny_beta
    skinny_aug_y_samples[t,] <- skinny_aug_y
    skinny_z_samples[t,] <- skinny_z
    skinny_w_samples[t,] <- skinny_w
  }
  return(list('exact_beta'=exact_beta_samples, 'exact_aug_y'=exact_aug_y_samples, 
              'exact_z'=exact_z_samples, 'exact_w'=exact_w_samples,
              'skinny_beta'=skinny_beta_samples, 'skinny_aug_y'=skinny_aug_y_samples, 
              'skinny_z'=skinny_z_samples, 'skinny_w'=skinny_w_samples))
}


