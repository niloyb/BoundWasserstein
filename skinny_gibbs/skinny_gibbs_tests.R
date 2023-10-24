rm(list = ls())
source(file = 'skinny_gibbs/skinny_gibbs_functions.R')

# Testing the beta update function
iterations <- 1000
test_beta_update_1 <- rep(NA,iterations)
for(i in 1:iterations){
  n <- 1 + round(exp(1))
  p <- 1 + round(2*exp(1))
  X <- matrix(runif(n*p),n,p)
  aug_y <- rnorm(n)
  w <- rexp(n)
  z <- rbinom(p,size = 1,prob = 0.5)
  tau0 <- runif(1)
  tau1 <- runif(1)
  
  beta_mean <- update_beta_exact(w, z, aug_y, X, tau0=tau0, tau1=tau1, u=rep(0,p), delta = rep(0,n))
  wX <- X*(w^0.5)
  wX_transpose <- t(wX)
  w_aug_y <- aug_y*(w^0.5)
  d <- z/(tau1^2)+(1-z)/(tau0^2)
  
  test_beta_update_1[i] <- sum((solve(t(wX)%*%wX+diag(d))%*%t(wX)%*%w_aug_y-beta_mean)^2)^0.5
}

max(test_beta_update_1)

# Testing the aug_y update function
iterations <- 5
test_aug_y_update_1 <- rep(NA,iterations)
test_aug_y_update_2 <- rep(NA,iterations)
for(i in 1:iterations){
  n <- 1
  p <- 1 + round(2*exp(1))
  X <- matrix(runif(n*p),n,p)
  aug_y <- rnorm(n)
  w <- rexp(n)
  beta <- rnorm(p)
  e <- 0
  
  sample1 <- sapply(c(1:100000), function(x){TruncatedNormal::rtnorm(1, mu=cpp_mat_vec_prod(X, beta), sd=1/sqrt(w), lb=-Inf, ub=0)})
  sample2 <- sapply(c(1:100000), function(x){update_aug_y(beta, w, e, X)})
  mean1 <- mean(sample1)
  mean2 <- mean(sample2)
  
  var1 <- var(sample1)
  var2 <- var(sample2)
  
  test_aug_y_update_1[i] <- abs(mean1-mean2)
  test_aug_y_update_2[i] <- abs(var1-var2)
}
max(test_aug_y_update_1)
max(test_aug_y_update_2)


# Testing the z update function
# z update
n <- 50
p <- 20
s <- 5
true_beta <- matrix(0,p,1)
true_beta[1:s] = 2^(-(seq(s)/4-9/4))
w <- rep(1,n)
tau0 <- 0.01
tau1 <- 1
q <- 0.05

update_z_exact(true_beta+0.01*rnorm(p), w, tau0, tau1, q)




n <- 10 #1 + round(exp(1))
p <- 20 #n + round(2*exp(1))
X <- matrix(runif(n*p),n,p)
aug_y <- rnorm(n)
w <- rexp(n)
z <- rbinom(p,size = 1,prob = 0.5)
tau0 <- runif(1)
tau1 <- runif(1)

# Testing the skinny beta update function
iterations <- 1000
test_beta_update_2 <- rep(NA,iterations)
for(i in 1:iterations){
  n <- 1 + round(exp(1))
  p <- 1 + round(2*exp(1))
  X <- matrix(runif(n*p),n,p)
  X <- matrix(scale(X),n,p)
  aug_y <- rnorm(n)
  w <- rexp(n)
  z <- rbinom(p,size = 1,prob = 0.5)
  z[1] <- 1
  tau0 <- runif(1)
  tau1 <- runif(1)
  
  beta_mean <- update_beta_skinny(w, z, aug_y, X, tau0=tau0, tau1=tau1, u=rep(0,p), delta = rep(0,n))
  wX <- X*(w^0.5)
  wX_transpose <- t(wX)
  w_aug_y <- aug_y*(w^0.5)
  d <- z/(tau1^2)+(1-z)/(tau0^2)
  
  diags <- diag(c(d[z==1]),nrow = sum(z==1),ncol = sum(z==1))
  
  test_beta_update_2[i] <- 
    sum((solve(t(wX[,z==1,drop=FALSE])%*%wX[,z==1,drop=FALSE]+diags)%*%
           t(wX[,z==1,drop=FALSE])%*%w_aug_y-beta_mean[z==1])^2)^0.5
}

max(test_beta_update_2)


# Testing the skinny z update function
n <- 100 #1 + round(exp(1))
p <- 200 #n + round(2*exp(1))
X <- matrix(runif(n*p),n,p)
X <- matrix(scale(X),n,p) 
aug_y <- rnorm(n)
s <- pi/sqrt(3) # sd of logistic distribution
w <- rep(1,n)/(s^2)
q <- 0.05
z <- rbinom(p,size = 1,prob = q)
tau0 <- 1/sqrt(n) # Following skinny gibbs paper
tau1 <- sqrt(max(1, p^(2.1)/(100*n))) # Following skinny gibbs paper
beta <- c(rnorm(10)*tau1, rnorm(p-10)*tau0)
X_scaled=TRUE

update_z_skinny <- function(beta, w, z, aug_y, X, tau0, tau1, q, 
                            u_crn=NULL, X_scaled=TRUE){
  p <- length(beta)
  if(is.null(u_crn)){u_crn <- runif(p)}
  
  log_prob1 <- log(q)+dnorm(beta, sd=tau1, log = TRUE)
  log_prob0 <- log(1-q)+dnorm(beta, sd=tau0, log=TRUE)
  
  active <- z==1
  X_active <- X[,active,drop=FALSE]
  beta_active <- beta[active]
  residuals <- aug_y-cpp_mat_vec_prod(X_active, beta_active)
  for (j in 1:p){
    
    Xjbetaj <- beta[j]*X[,j]
    residuals_j <- residuals + Xjbetaj
    if(X_scaled){
      modification <- # ONLY when w has one unique value and X is scaled
        sum(Xjbetaj*residuals_j)*w[1] + 0.5*beta[j]^2*(n-1)*(1-w[1])
    } else{
      modification <- 
        sum(Xjbetaj*w*residuals_j) + 0.5*beta[j]^2*sum(X[,j]*(1-w)*X[,j])
    }
    
    # active[j] <- FALSE
    # X_active_minus_j <- X[,active_minus_j,drop=FALSE]
    # beta_active_minus_j <- beta[active_minus_j]
    # residuals_j <- aug_y-cpp_mat_vec_prod(X_active_minus_j, beta_active_minus_j)
    # modification <- 
    #   sum((beta[j]*X[,j])*w*residuals_j) + 0.5*beta[j]^2*sum(X[,j]*(1-w)*X[,j])
    
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
colMeans(sapply(c(1:100),function(x){update_z_skinny(beta, w, z, aug_y, X, tau0, tau1, q)}))


# Update aug_y crn
n <- 50
p <- 200
s0 <- 10
true_beta <- matrix(0,p,1)
true_beta[1:s0] = 2^(-(seq(s0)/4-9/4))
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
true_aug_y = rlogis(n, location = X%*%true_beta) # variance pi^2/3
y <- ifelse(true_aug_y>0,1,0) # Logistic response
X <- matrix(scale(X),n,p) ### important for Bayesian variable selection
s <- pi/sqrt(3) # sd of logistic distribution
w <- rep(1,n)/(s^2)
beta <- rnorm(p)

update_aug_y(beta, w, y, X)



# Update z crn
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
# Choice of q, tau0, tau1
K <- max(10,log(n))
q_seq <- seq(0,1,0.0001)
probs <- abs(pbinom(K,p,q_seq)-0.9)
q_index <- which(probs==min(probs))
q <- q_seq[q_index] # Following skinny gibbs paper
tau0 <- 1/sqrt(n) # Following skinny gibbs paper
tau1 <- sqrt(max(1, p^(2.1)/(100*n))) # Following skinny gibbs paper
update_z_exact(beta, w, tau0, tau1, q)




# Coupled kernel test
n <- 50
p <- 200
s0 <- 10
true_beta <- matrix(0,p,1)
true_beta[1:s0] = 2^(-(seq(s0)/4-9/4))
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
true_aug_y = rlogis(n, location = X%*%true_beta) # variance pi^2/3
y <- ifelse(true_aug_y>0,1,0) # Logistic response
X <- matrix(scale(X),n,p) ### important for Bayesian variable selection

exact_beta <- rnorm(p) 
exact_aug_y <- rnorm(n)
exact_z <- rbinom(p,1,0.5)
s <- pi/sqrt(3) # sd of logistic distribution
exact_w <- rep(1,n)/(s^2)
skinny_beta <- rnorm(p)
skinny_aug_y <- rnorm(n)
skinny_z <- rbinom(p,1,0.5)
skinny_w <- rep(1,n)/(s^2)

# Choice of q, tau0, tau1
K <- max(10,log(n))
q_seq <- seq(0,1,0.0001)
probs <- abs(pbinom(K,p,q_seq)-0.9)
q_index <- which(probs==min(probs))
q <- q_seq[q_index] # Following skinny gibbs paper
tau0 <- 1/sqrt(n) # Following skinny gibbs paper
tau1 <- sqrt(max(1, p^(2.1)/(100*n))) # Following skinny gibbs paper
X_scaled=TRUE

exact_skinny_coupled_kernel(exact_beta, exact_aug_y, exact_z, exact_w, 
                            skinny_beta, skinny_aug_y, skinny_z, skinny_w, 
                            X, y, tau0, tau1, q, X_scaled)


#### Coupled chain test
n <- 100
p <- 1000
s0 <- 10
true_beta <- matrix(0,p,1)
true_beta[1:s0] = 2^(-(seq(s0)/4-9/4))
X <- matrix(rnorm(n*p), nrow = n, ncol = p)
# Error terms
# error_std <- 0.1
# error_terms = error_std*rnorm(n, mean = 0, sd = 1)
# true_aug_y = X%*%true_beta + error_terms
true_aug_y = rlogis(n, location = X%*%true_beta) # variance pi^2/3
y <- ifelse(true_aug_y>0,1,0) # Logistic response
X <- matrix(scale(X),n,p) ### important for Bayesian variable selection

# Choice of q, tau0, tau1
K <- max(10,log(n))
q_seq <- seq(0,1,0.0001)
probs <- abs(pbinom(K,p,q_seq)-0.9)
q_index <- which(probs==min(probs))
q <- q_seq[q_index] # Following skinny gibbs paper
tau0 <- 1/sqrt(n) # Following skinny gibbs paper
tau1 <- sqrt(max(1, p^(2.1)/(100*n))) # Following skinny gibbs paper

# Fixing w
s <- pi/sqrt(3) # sd of logistic distribution
w <- rep(1,n)/(s^2)

burnin <- 5e1
chain_length <- 1e2

# Testing skinny z update
coupled_chain_test <- exact_skinny_chain(X, y, tau0, tau1, q, w, chain_length=1e3, rinit=NULL,X_scaled=TRUE)

