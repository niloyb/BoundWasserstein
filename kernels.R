#### #### #### #### #### #### #### 
#### ULA and MALA Functions #### 
#### #### #### #### #### #### 
# Single chain kernels
ula_proposal_mean <- function(x, grad_log_pdf, sigma_mh){
  return(x + 0.5*(sigma_mh^2)*grad_log_pdf(x))
}
# ULA
kernel_ula <- function(current_state, grad_log_pdf,sigma_mh){
  rw_step <- rnorm(length(current_state))
  current_state <- 
    ula_proposal_mean(current_state, grad_log_pdf, sigma_mh) + sigma_mh*rw_step
  # current_state <-   current_state + 
  #   0.5*(sigma_mh^2)*grad_log_pdf(current_state) + sigma_mh*rw_step
  return(current_state)
}
ula <- function(init_state, grad_log_pdf, sigma_mh, iterations){
  if (length(sigma_mh)==1){sigma_mh <- rep(sigma_mh, iterations)}
  samples <- matrix(NA, nrow = iterations, ncol = length(init_state))
  for(i in 1:iterations){
    samples[i,] <- init_state
    init_state <- kernel_ula(init_state, grad_log_pdf, sigma_mh[i])
  }
  return(samples)
}

# MALA
kernel_mala <- function(current_state, log_pdf, grad_log_pdf, sigma_mh){
  rw_step <- rnorm(length(current_state))
  proposed_state <- 
    ula_proposal_mean(current_state, grad_log_pdf, sigma_mh) + sigma_mh*rw_step
  # proposed_state <- current_state + 0.5*(sigma_mh^2)*grad_log_pdf(current_state) + sigma_mh*rw_step
  current_state <- as.vector(current_state)
  proposed_state <- as.vector(proposed_state)
  # accept <- log(runif(1)) < min(0, log_pdf(proposed_state)-log_pdf(current_state) +
  #                               mvtnorm::dmvnorm(current_state, mean = proposed_state + 0.5*(sigma_mh^2)*grad_log_pdf(proposed_state), sigma = sigma_mh^2*diag(length(current_state)), log = TRUE) -
  #                               mvtnorm::dmvnorm(proposed_state, mean = current_state + 0.5*(sigma_mh^2)*grad_log_pdf(current_state), sigma = sigma_mh^2*diag(length(current_state)), log = TRUE))
  accept <- log(runif(1)) < log_pdf(proposed_state)-log_pdf(current_state) +
    -0.5*sum((current_state-ula_proposal_mean(proposed_state, grad_log_pdf, sigma_mh))^2)/sigma_mh^2+
    +0.5*sum((proposed_state-ula_proposal_mean(current_state, grad_log_pdf, sigma_mh))^2)/sigma_mh^2
  if(accept){current_state <- proposed_state}
  return(current_state)
}
mala <- function(init_state, log_pdf, grad_log_pdf, sigma_mh, iterations){
  if (length(sigma_mh)==1){sigma_mh <- rep(sigma_mh, iterations)}
  samples <- matrix(NA, nrow = iterations, ncol = length(init_state))
  for(i in 1:iterations){
    samples[i,] <- init_state
    init_state <- kernel_mala(init_state, log_pdf, grad_log_pdf, sigma_mh[i])
  }
  return(samples)
}

# Coupled Kernels
# Reflection Coupling of N(mu1, Id) and N(mu2, Id)
reflection_couping_normal_Idcov <- function(mu1, mu2, norm_crn=NULL, unif_crn=NULL){
  if(is.null(norm_crn)){
    Xdot <- rnorm(length(mu1))
  } else {
    Xdot <- norm_crn
  }
  if(is.null(unif_crn)){
    uniform_crn <- runif(1)
  } else {
    uniform_crn <- unif_crn
  }
  z <- mu1-mu2
  if(sum(z^2)==0){return(cbind(mu1+Xdot, mu1+Xdot))}
  e <- z/sum(z^2)^0.5
  if(log(uniform_crn) < -0.5*(sum((Xdot+z)^2)-sum((Xdot)^2))){
    Ydot <- Xdot + z
  } else{
    Ydot <- Xdot - 2*sum(e*Xdot)*e
  }
  return(cbind((mu1+Xdot), (mu2+Ydot)))
}
# Reflection Coupling of N(mu1, Sigma1) and N(mu2, Sigma2)
reflection_couping_normal <- 
  function(mu1, mu2, cov_mat1_cholesky, cov_mat2_cholesky,
           cov_mat1_cholesky_inv, cov_mat2_cholesky_inv,
           norm_crn=NULL, unif_crn=NULL){
    mu1tilde <- cov_mat1_cholesky_inv%*%mu1
    mu2tilde <- cov_mat2_cholesky_inv%*%mu2
    residuals <- 
      reflection_couping_normal_Idcov(mu1tilde, mu2tilde, norm_crn, unif_crn)
    return(cbind((mu1+cov_mat1_cholesky%*%residuals[,1]), 
                 (mu2+cov_mat2_cholesky%*%residuals[,2])))
  }

# Coupled ULA
coupled_kernel_ula <-
  function(current_pi, current_q, pi_grad_log_pdf, q_grad_log_pdf,
           sigma_mh, reflect_threshold=0, norm_crn=NULL, unif_crn=NULL){
    if(is.null(norm_crn)){norm_crn <- rnorm(length(current_q))}
    if(is.null(unif_crn)){unif_crn <- runif(1)}
    
    current_q <- as.vector(current_q)
    current_pi <- as.vector(current_pi)
    r <- sum((current_q-current_pi)^2)^0.5
    if (length(sigma_mh)==1){
      sigma_mh_q <- sigma_mh
      sigma_mh_pi <- sigma_mh
    } else{
      sigma_mh_q <- sigma_mh[1]
      sigma_mh_pi <- sigma_mh[2]
    }
    
    if (r < reflect_threshold){
      # Reflection coupling when close
      mu_q <- (current_q + 0.5*(sigma_mh_q^2)*q_grad_log_pdf(current_q))/sigma_mh_q
      mu_pi <- (current_pi + 0.5*(sigma_mh_pi^2)*pi_grad_log_pdf(current_pi))/sigma_mh_pi
      coupled_output <- reflection_couping_normal_Idcov(mu_q, mu_pi, norm_crn, unif_crn)
      current_q <- coupled_output[,1]*sigma_mh_q
      current_pi <- coupled_output[,2]*sigma_mh_pi
      return(list('P'=current_pi, 'Q'=current_q))
    } else{
      # Synchronous coupling when far away
      rw_step <- norm_crn
      current_q <- current_q + 0.5*(sigma_mh_q^2)*q_grad_log_pdf(current_q) +
        rw_step*sigma_mh_q
      current_pi <- current_pi + 0.5*(sigma_mh_pi^2)*pi_grad_log_pdf(current_pi) +
        rw_step*sigma_mh_pi
      return(list('P'=current_pi, 'Q'=current_q))
    }
  }
coupled_ula <- function(init_pi, init_q, pi_grad_log_pdf, q_grad_log_pdf,
                        sigma_mh, iterations, reflect_threshold=0){
  if (length(sigma_mh)==1){sigma_mh_mat <- as.matrix(rep(sigma_mh, iterations))}
  if (length(sigma_mh)==2){
    sigma_mh_mat <- matrix(NA, nrow = iterations, ncol = 2)
    sigma_mh_mat[,1] <- sigma_mh[1]
    sigma_mh_mat[,2] <- sigma_mh[2]
  }
  samples_q <- matrix(NA, nrow = iterations, ncol = length(init_q))
  samples_pi <- matrix(NA, nrow = iterations, ncol = length(init_pi))
  for(i in 1:iterations){
    samples_q[i,] <- init_q
    samples_pi[i,] <- init_pi
    output <- 
      coupled_kernel_ula(init_pi, init_q, pi_grad_log_pdf, q_grad_log_pdf,
                         sigma_mh_mat[i,], reflect_threshold)
    init_pi <- output$P
    init_q <- output$Q
  }
  return(list('P'=samples_pi, 'Q'=samples_q))
}

# Coupled MALA kernel
coupled_kernel_mala <- function(current_pi, current_q, pi_log_pdf, q_log_pdf,
                                pi_grad_log_pdf, q_grad_log_pdf, sigma_mh, 
                                reflect_threshold=0, pi_correction=TRUE,
                                q_correction=TRUE,
                                norm_crn=NULL, unif_crn1=NULL, unif_crn2=NULL){
  current_q <- as.vector(current_q)
  current_pi <- as.vector(current_pi)
  proposals <- 
    coupled_kernel_ula(current_pi, current_q, pi_grad_log_pdf, q_grad_log_pdf,
                       sigma_mh, reflect_threshold, norm_crn, unif_crn1)
  proposal_pi <- proposals$P
  proposal_q <- proposals$Q
  if(is.null(unif_crn2)){unif_crn2 <- runif(1)}
  log_uniform <- log(unif_crn2)
  
  if (length(sigma_mh)==1){
    sigma_mh_q <- sigma_mh
    sigma_mh_pi <- sigma_mh
  } else{
    sigma_mh_q <- sigma_mh[1]
    sigma_mh_pi <- sigma_mh[2]
  }
  
  if(q_correction){
    # accept_q <- 
    #   log_uniform < min(0, q_log_pdf(proposal_q)-q_log_pdf(current_q) +
    #                       mvtnorm::dmvnorm(current_q, mean = proposal_q + 0.5*(sigma_mh_q^2)*q_grad_log_pdf(proposal_q), sigma = sigma_mh_q^2*diag(length(current_q)), log = TRUE) -
    #                       mvtnorm::dmvnorm(proposal_q, mean = current_q + 0.5*(sigma_mh_q^2)*q_grad_log_pdf(current_q), sigma = sigma_mh_q^2*diag(length(current_q)), log = TRUE))
    accept_q <- log_uniform < q_log_pdf(proposal_q)-q_log_pdf(current_q) +
      -0.5*sum((current_q-ula_proposal_mean(proposal_q, q_grad_log_pdf, sigma_mh_q))^2)/sigma_mh_q^2+
      +0.5*sum((proposal_q-ula_proposal_mean(current_q, q_grad_log_pdf, sigma_mh_q))^2)/sigma_mh_q^2
  } else {
    accept_q <- TRUE
  }
  if(pi_correction){
    # accept_pi <- 
    #   log_uniform < min(0, pi_log_pdf(proposal_pi)-pi_log_pdf(current_pi) +
    #                       mvtnorm::dmvnorm(current_pi, mean = proposal_pi + 0.5*(sigma_mh_pi^2)*pi_grad_log_pdf(proposal_pi), sigma = sigma_mh_pi^2*diag(length(current_pi)), log = TRUE) -
    #                       mvtnorm::dmvnorm(proposal_pi, mean = current_pi + 0.5*(sigma_mh_pi^2)*pi_grad_log_pdf(current_pi), sigma = sigma_mh_pi^2*diag(length(current_pi)), log = TRUE))
    accept_pi <- log_uniform < pi_log_pdf(proposal_pi)-pi_log_pdf(current_pi) +
      -0.5*sum((current_pi-ula_proposal_mean(proposal_pi, pi_grad_log_pdf, sigma_mh_pi))^2)/sigma_mh_pi^2+
      +0.5*sum((proposal_pi-ula_proposal_mean(current_pi, pi_grad_log_pdf, sigma_mh_pi))^2)/sigma_mh_pi^2
  } else {
    accept_pi <- TRUE
  }
  if(accept_q){current_q <- proposal_q}
  if(accept_pi){current_pi <- proposal_pi}
  return(list('P'=current_pi, 'Q'=current_q))
}
coupled_mala <- function(init_pi, init_q, pi_log_pdf, q_log_pdf,
                         pi_grad_log_pdf, q_grad_log_pdf, sigma_mh,
                         iterations, reflect_threshold=0,
                         pi_correction=TRUE, q_correction=TRUE){
  if (length(sigma_mh)==1){sigma_mh_mat <- as.matrix(rep(sigma_mh, iterations))}
  if (length(sigma_mh)==2){
    sigma_mh_mat <- matrix(NA, nrow = iterations, ncol = 2)
    sigma_mh_mat[,1] <- sigma_mh[1]
    sigma_mh_mat[,2] <- sigma_mh[2]
  }
  samples_q <- matrix(NA, nrow = iterations, ncol = length(init_q))
  samples_pi <- matrix(NA, nrow = iterations, ncol = length(init_pi))
  for(i in 1:iterations){
    samples_q[i,] <- init_q
    samples_pi[i,] <- init_pi
    output <- 
      coupled_kernel_mala(init_pi, init_q, pi_log_pdf, q_log_pdf,
                          pi_grad_log_pdf, q_grad_log_pdf, sigma_mh_mat[i,],
                          reflect_threshold, pi_correction, q_correction)
    init_pi <- output$P
    init_q <- output$Q
  }
  return(list('Q'=samples_q, 'P'=samples_pi))
}

#### #### #### #### ####
#### HMC Functions #### 
#### #### #### #### ####

# Leapfrog integrator function
leapfrog <- function(x, v, gradlogtarget, stepsize, nsteps) {
  v <- v + stepsize * gradlogtarget(x)/2
  for (step in 1:nsteps) {
    x <- x + stepsize * v
    if (step != nsteps) {
      v <- v + stepsize * gradlogtarget(x)
    }
  }
  v <- v + stepsize * gradlogtarget(x)/2
  return(list(x = x, v = v))
}
# Single chain HMC kernel
hmc_kernel <- function(current_state, current_pdf, logtarget, 
                       gradlogtarget, stepsize, nsteps) {
  current_v <- rnorm(dimension)
  leapfrog_result <- 
    leapfrog(current_state, current_v, gradlogtarget, stepsize, nsteps)
  proposed_v <- -leapfrog_result$v
  proposed_x <- leapfrog_result$x
  proposed_pdf <- logtarget(proposed_x)
  accept_ratio <- proposed_pdf - current_pdf
  accept_ratio <- accept_ratio + sum(current_v^2)/2 - sum(proposed_v^2)/2
  accept <- FALSE
  if (is.finite(accept_ratio)) {
    accept <- (log(runif(1)) < accept_ratio)
  }
  if (accept) {
    current_state <- proposed_x
    current_pdf <- proposed_pdf
    accept <- TRUE
  }
  return(list(chain_state = current_state, chain_pdf = current_pdf, 
              accept = accept))
}
# Single chain HMC
hmc <- function(init_state, init_pdf, logtarget, gradlogtarget, 
                stepsize, nsteps, chain_length){
  hmc_chain <- matrix(nrow = chain_length, ncol = length(init_state))
  current_state <- init_state
  current_pdf <- init_pdf
  for (t in 1:chain_length){
    hmc_step <- hmc_kernel(current_state, current_pdf, logtarget, gradlogtarget, 
                           stepsize, nsteps)
    current_state <- hmc_step$chain_state
    current_pdf <- hmc_step$chain_pdf
    hmc_chain[t,] <- current_state
  }
  return(list(chain=hmc_chain))
}

# Coupled HMC kernel
coupled_hmc_kernel <- 
  function(current_pi, current_q, current_pdf_pi, current_pdf_q, 
           pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf, 
           pi_stepsize, q_stepsize, pi_nsteps, q_nsteps, reflect_threshold,
           pi_correction, q_correction) {
  r <- sum((current_pi-current_q)^2)^0.5
  if (r < reflect_threshold){
    # Reflection Coupling of momentum when close from current_pi, current_q
    coupled_v <- reflection_couping_normal_Idcov(current_pi, current_q)
    coupled_v <- coupled_v - cbind(current_pi, current_q)
  } else {
    # Synchronous coupling of momentum when far away
    coupled_v <- rnorm(length(current_pi))
    coupled_v <- cbind(coupled_v,coupled_v)
  }
  
  leapfrog_result1 <- 
    leapfrog(current_pi, coupled_v[,1], pi_grad_log_pdf, pi_stepsize, pi_nsteps)
  leapfrog_result2 <- 
    leapfrog(current_q, coupled_v[,2], q_grad_log_pdf, q_stepsize, q_nsteps)
  proposed_v1 <- -leapfrog_result1$v
  proposed_x1 <- leapfrog_result1$x
  proposed_pdf1 <- pi_log_pdf(proposed_x1)
  proposed_v2 <- -leapfrog_result2$v
  proposed_x2 <- leapfrog_result2$x
  proposed_pdf2 <- q_log_pdf(proposed_x2)
  
  log_unif_crn <- log(runif(1))
  if(pi_correction){ # Metropolis Correction
    accept_ratio1 <- proposed_pdf1 - current_pdf_pi
    kinetic_energy1 <- sum(coupled_v[,1]^2)/2
    accept_ratio1 <- accept_ratio1 + kinetic_energy1 - sum(proposed_v1^2)/2
    if(is.finite(accept_ratio1)&(log_unif_crn < accept_ratio1)){
      current_pi <- proposed_x1
      current_pdf_pi <- proposed_pdf1
    }
  } else { # No Metropolis Correction
    current_pi <- proposed_x1
    current_pdf_pi <- proposed_pdf1
  }
  if(q_correction){ # Metropolis Correction
    accept_ratio2 <- proposed_pdf2 - current_pdf_pi
    kinetic_energy2 <- sum(coupled_v[,2]^2)/2
    accept_ratio2 <- accept_ratio2 + kinetic_energy2 - sum(proposed_v2^2)/2
    if(is.finite(accept_ratio2)&(log_unif_crn < accept_ratio2)){
      current_q <- proposed_x2
      current_pdf_q <- proposed_pdf2
    }
  } else { # No Metropolis Correction
    current_q <- proposed_x2
    current_pdf_q <- proposed_pdf2
  }
  return(list(current_pi = current_pi, current_q = current_q, 
              current_pdf_pi = current_pdf_pi, current_pdf_q = current_pdf_q))
}

# HMC with synchronous coupling
coupled_hmc <- 
  function(init_pi, init_q, init_pdf_pi, init_pdf_q, 
           pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf,
           pi_stepsize, q_stepsize, pi_nsteps, q_nsteps,
           chain_length, reflect_threshold=0,
           pi_correction=TRUE, q_correction=TRUE){
  hmc_chain1 <- matrix(nrow = chain_length, ncol = length(init_pi))
  hmc_chain2 <- matrix(nrow = chain_length, ncol = length(init_q))
  current_pi <- init_pi
  current_q <- init_q
  current_pdf_pi <- init_pdf_pi
  current_pdf_q <- init_pdf_q
  hmc_chain1[1,] <- current_pi
  hmc_chain2[1,] <- current_q
  for (t in 1:(chain_length-1)){
    coupled_hmc_step <- 
      coupled_hmc_kernel(current_pi, current_q, current_pdf_pi, current_pdf_q, 
                         pi_log_pdf, q_log_pdf, pi_grad_log_pdf, q_grad_log_pdf, 
                         pi_stepsize, q_stepsize, pi_nsteps, q_nsteps, 
                         reflect_threshold, pi_correction, q_correction)
    current_pi <- coupled_hmc_step$current_pi
    current_q <- coupled_hmc_step$current_q
    current_pdf_pi <- coupled_hmc_step$current_pdf_pi
    current_pdf_q <- coupled_hmc_step$current_pdf_q
    hmc_chain1[(t+1),] <- current_pi
    hmc_chain2[(t+1),] <- current_q
  }
  return(list('P'=hmc_chain1, 'Q'=hmc_chain2))
}






