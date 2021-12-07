# Estimators
require(doParallel)
require(Rcpp)
require(RcppEigen)

#### Wp upper bounds ####
wp_ub_estimate <- 
  function(coupled_chain_sampler, no_chains, p=1, metric, parallel=TRUE){
    if(parallel){
      trajectories <- foreach(i = c(1:no_chains), .combine = rbind)%dopar%{
        chain <- coupled_chain_sampler()
        return(metric(chain$Q, chain$P))
      }
    } else {
      trajectories <- foreach(i = c(1:no_chains), .combine = rbind)%do%{
        chain <- coupled_chain_sampler()
        print(i)
        return(metric(chain$Q, chain$P))
      }
    }
    
    if(no_chains==1){
      wp_power_p_ub_mean_tracjectory <- trajectories^p
      wp_ub_se <- NA
    } else{
      wp_power_p_ub_mean_tracjectory <- colMeans(trajectories^p)
      
      wp_power_p_ubs <- rowMeans(trajectories^p)
      wp_ub_se <- (sd(wp_power_p_ubs)/sqrt(no_chains))^(1/p)
      
      # if(p==1){
      #   wp_ub_se <- sd(rowMeans(trajectories))
      # } else{
      #   wpub_fc <- function(trajectories, i){return(mean(trajectories[i,]^p)^(1/p))}
      #   boot_ub <- boot::boot(trajectories, wpub_fc, R=boot_sample)
      #   wp_ub_se <- sd(boot_ub$t)
      # }
    }
    wp_power_p_ub <- mean(wp_power_p_ub_mean_tracjectory)
    
    return(list('wp_power_p_ub_tracjectories' = trajectories^p,
                'wp_power_p_ub_mean_tracjectory' = wp_power_p_ub_mean_tracjectory,
                'wp_power_p_ub'=wp_power_p_ub, 'wp_ub'=wp_power_p_ub^(1/p),
                'wp_ub_se'=wp_ub_se))
  }

#### 2-Wasserstein with L2 distance lower bound ####
## Rcpp functions ##
Rcpp::cppFunction('
Eigen::MatrixXd cpp_crossprod(const Eigen::MatrixXd X){
Eigen::MatrixXd output;
output = (X.transpose())*X;
return output;
}', depends = 'RcppEigen')

w2l2_lb_gelbrich <- function(chain1,chain2){ # LB of gelbrich
  no_data_points1 <- dim(chain1)[1]
  no_data_points2 <- dim(chain2)[1]
  mu1 <- colMeans(chain1)
  mu2 <- colMeans(chain2)
  chain1 <- sweep(chain1, 2, mu1)/sqrt(no_data_points1-1)
  chain2 <- sweep(chain2, 2, mu2)/sqrt(no_data_points2-1)
  
  # sqrtSigma1 equals sqrt root of covariance matrix t(chain1)%*%(chain1)
  chain1_svd <- svd(t(chain1))
  #sqrtSigma1 <- (chain1_svd$u)%*%diag((chain1_svd$d))%*%t(chain1_svd$u)
  sqrtSigma1 <- cpp_crossprod(t(chain1_svd$u)*((chain1_svd$d)^0.5)) # faster
  trace_Sigma1 <- sum((chain1_svd$d)^2)
  
  # sqrtSigma2 equals sqrt root of covariance matrix chain2%*%t(chain2)
  chain2_svd <- svd(t(chain2))
  # sqrtSigma2 <- (chain2_svd$u)%*%diag((chain2_svd$d))%*%t(chain2_svd$u)
  sqrtSigma2 <- cpp_crossprod(t(chain2_svd$u)*((chain2_svd$d)^0.5)) # faster
  trace_Sigma2 <- sum((chain2_svd$d)^2)
  
  # sqrtCrossSigma equals sqrt root of sqrtSigma1%*%Sigma2%*%sqrtSigma1
  # sqrtCrossSigma <- (CrossSigma_svd$u)%*%diag((CrossSigma_svd$d))%*%t(CrossSigma_svd$u)
  # sqrtCrossSigma <- cpp_crossprod(t(CrossSigma_svd$u)*((CrossSigma_svd$d)^0.5)) # faster
  sqrt_prod <- cpp_prod(sqrtSigma1, sqrtSigma2)
  CrossSigma_svd <- svd(sqrt_prod)
  trace_sqrtCrossSigma <- sum((CrossSigma_svd$d))
  
  out <- trace_Sigma1+trace_Sigma2-2*trace_sqrtCrossSigma + sum((mu1-mu2)^2)
  if(abs(out)<1e-10){out <- abs(out)} # If small output 0
  return(out^0.5)
}

## W2L2LB based on marginals ##
W2LBmarginals <- function(chain1, chain2){
  no_data_points1 <- dim(chain1)[1]
  no_data_points2 <- dim(chain2)[1]
  if(no_data_points1!=no_data_points2){stop("samples sizes must be equal")}
  
  dimension <- dim(chain1)[2]
  
  W2marginals <- rep(NA, dimension)
  for(i in 1:dimension){
    W2marginals[i] <- 
      mean((chain1[,i][order(chain1[,i])]-chain2[,i][order(chain2[,i])])^2)
  }
  return(sum(W2marginals)^0.5)
}


#### 2-Wasserstein with L2 distance upper and lower bounds ####
# L2 metric
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}
W2L2_UBLB_estimates <- 
  function(chain_sampler, no_chains, parallel=TRUE, lb='max_gelbrich_marginals'){
    if(parallel){
      output <- foreach(i = c(1:no_chains), .combine = rbind)%dopar%{
        chain <- chain_sampler()
        metric_l2_traj <- metric_l2(chain$Q, chain$P)
        lb_gelbrich <- NA
        if(grepl('gelbrich', lb, fixed = TRUE)){
          lb_gelbrich <- w2l2_lb_gelbrich(chain$Q, chain$P)
        }
        lb_marginals <- NA
        if(grepl('marginals', lb, fixed = TRUE)){
          lb_marginals <- W2LBmarginals(chain$Q, chain$P)
        }
        return(list(lb_gelbrich=lb_gelbrich, lb_marginals=lb_marginals,
                    metric_l2_traj=metric_l2_traj))
      }
    } else {
      output <- foreach(i = c(1:no_chains), .combine = rbind)%do%{
        chain <- chain_sampler()
        metric_l2_traj <- metric_l2(chain$Q, chain$P)
        lb_gelbrich <- NA
        if(grepl('gelbrich', lb, fixed = TRUE)){
          lb_gelbrich <- w2l2_lb_gelbrich(chain$Q, chain$P)
        }
        lb_marginals <- NA
        if(grepl('marginals', lb, fixed = TRUE)){
          lb_marginals <- W2LBmarginals(chain$Q, chain$P)
        }
        return(list(lb_gelbrich=lb_gelbrich, lb_marginals=lb_marginals,
                    metric_l2_traj=metric_l2_traj))
      }
    }
    
    if(no_chains==1){
      trajectories <- output$metric_l2_traj
      W2L2_UBs <- mean(trajectories^2)^0.5
      W2L2_UBmean <- W2L2_UBs
      W2L2_UBsd <- NA
      
      W2L2_LBs <- max(output$lb_gelbrich, output$lb_marginals, na.rm = TRUE)
      W2L2_LBmean <- W2L2_LBs
      W2L2_LBsd <- NA
    } else {
      trajectories <- do.call(rbind, output[,"metric_l2_traj"])
      W2L2_squared_UBs <- rowMeans(trajectories^2)
      W2L2_UBs <- W2L2_squared_UBs^0.5
      W2L2_UBmean <- mean(W2L2_squared_UBs)^0.5
      W2L2_UBsd <- (sd(W2L2_squared_UBs)/sqrt(no_chains))^0.5
      
      W2L2_squared_LBs <- 
        pmax(do.call(rbind, output[,"lb_gelbrich"]),
             do.call(rbind, output[,"lb_marginals"]), na.rm = TRUE)^2
      W2L2_LBs <- W2L2_squared_LBs^0.5
      W2L2_LBmean <- mean(W2L2_squared_LBs)^0.5
      W2L2_LBsd <- (sd(W2L2_squared_LBs)/sqrt(no_chains))^0.5
    }

    return(list(W2L2_UBs=W2L2_UBs,'W2L2_UBmean'=W2L2_UBmean,'W2L2_UBsd'=W2L2_UBsd,
                W2L2_LBs=W2L2_LBs,'W2L2_LBmean'=W2L2_LBmean,'W2L2_LBsd'=W2L2_LBsd))
  }



