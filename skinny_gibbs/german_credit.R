rm(list = ls())
setwd('/Users/niloybiswas/Dropbox/Apps/Overleaf/couplings/code/paper_examples/')
source(file = 'skinny_gibbs/skinny_gibbs_functions.R')
source(file = 'estimators.R')


data <- read.csv('skinny_gibbs/german_credit.csv')
y <- data[,"Creditability"]

x.categorical <- c('Account.Balance', 'Payment.Status.of.Previous.Credit', 'Purpose', 'Value.Savings.Stocks',
                   'Length.of.current.employment', 'Sex...Marital.Status', 'Guarantors', 'Most.valuable.available.asset',
                   'Concurrent.Credits', 'Type.of.apartment', 'Occupation', 'Telephone', 'Foreign.Worker')
x.quant <- c('Duration.of.Credit..month.', 'Credit.Amount', 'Instalment.per.cent', 'Duration.in.Current.address',
             'Age..years.', 'No.of.Credits.at.this.Bank', 'No.of.dependents')
for(x in x.categorical){
  data[,x] = as.factor(data[,x])
}
fmla <- paste('~',paste(c(x.quant,x.categorical),collapse ='+'))
X <- model.matrix(formula(fmla), data=data)

X <- X[,-1] # remove intercept 
n <- nrow(X)
p <- ncol(X)
X <- matrix(scale(X),n,p) 

# Choice of q, tau0, tau1: following skinny Gibbs paper
K <- max(10,log(n))
q_seq <- seq(0,1,0.0001)
probs <- abs(pbinom(K,p,q_seq)-0.9)
q_index <- which(probs==min(probs))
q <- q_seq[q_index]
tau0 <- 1/sqrt(n)
tau1 <- sqrt(max(1, p^(2.1)/(100*n)))

# Fixing w
s <- pi/sqrt(3) # sd of logistic distribution
w <- rep(1,n)/(s^2)

burnin <- 1
chain_length <- burnin+5e2
#no_chains <- 20
# # Exact chain
exact_chain <- full_gibbs(X, y, tau0, tau1, q, w, chain_length=chain_length)
which(colMeans(exact_chain$z[(burnin:chain_length),])>0.5)
# matplot(exact_chain$beta[(burnin:chain_length),], type='l')






# Interactions
interaction_terms <- matrix(nrow = n, ncol = p*(p-1) / 2)
index <- 1
for (j in 1:(p-1)){
  for (jprime in (j+1):p){
    interaction_terms[, index] <- X[, j] * X[, jprime]
    index <- index + 1
  }
}

X_interaction <- cbind(X, scale(interaction_terms))
colnames(X_interaction) <- NULL
#dimension <- ncol(design_matrix) + 2
n <- dim(X_interaction)[1]
p <- dim(X_interaction)[2]
X_interaction <- matrix(scale(X_interaction),n,p)

# Choice of q, tau0, tau1: following skinny gibbs paper
K <- max(10,log(n))
q_seq <- seq(0,1,0.0001)
probs <- abs(pbinom(K,p,q_seq)-0.9)
q_index <- which(probs==min(probs))
q <- q_seq[q_index]
tau0 <- 1/sqrt(n)
tau1 <- sqrt(max(1, p^(2.1)/(100*n)))

# Fixing w
s <- pi/sqrt(3) # sd of logistic distribution
w <- rep(1,n)/(s^2)

burnin <- 1
chain_length <- burnin+1e1
# no_chains <- 20
# # Exact chain
exact_chain <- full_gibbs(X, y, tau0, tau1, q, w, chain_length=chain_length)
# # which(colMeans(exact_chain$z[(burnin:chain_length),])>0.5)
# matplot(exact_chain$beta[(burnin:chain_length),], type='l')
# # Skinny chain
skinny_chain <- skinny_gibbs(X_interaction, y, tau0, tau1, q, w, chain_length=chain_length)
# # which(colMeans(skinny_chain$z[(burnin:chain_length),])>0.5)
# matplot(skinny_chain$beta[(burnin:chain_length),], type='l')

# Metrics considered: L2 and Hamming. 
metric_l2 <- function(x,y){
  if(is.null(dim(x))){return(sum((x-y)^2)^0.5)} else {return(rowSums((x-y)^2)^0.5)}
}
# Note: L2 distance squared now gives the Hamming distance
metric_hamming <- function(x,y){metric_l2(x,y)^2}

crn_exact_chain <- coupled_exact_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
crn_skinny_chain <- coupled_skinny_chain(X, y, tau0, tau1, q, w, chain_length=chain_length)
matplot(metric_hamming(crn_exact_chain$z1,crn_exact_chain$z2), type = 'l')
matplot(metric_hamming(crn_skinny_chain$z1,crn_skinny_chain$z2), type = 'l')

matplot(metric_l2(crn_exact_chain$beta1,crn_exact_chain$beta2), type = 'l')
matplot(metric_l2(crn_skinny_chain$beta1,crn_skinny_chain$beta2), type = 'l')



matplot(crn_exact_chain$beta1[(burnin:chain_length),], type='l')
matplot(crn_exact_chain$beta2[(burnin:chain_length),], type='l')
matplot(crn_skinny_chain$beta1, type='l')
matplot(crn_skinny_chain$beta2, type='l')

matplot(metric_hamming(crn_skinny_chain$beta1,crn_skinny_chain$beta2), type = 'l')


