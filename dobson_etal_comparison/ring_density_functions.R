

# Example 1
# PDF up to normalizing constant
V_ring <- function(x){
  return((x[1]^2+x[2]^2-1)^2)
}
ring_logpdf <- function(x, sigma){
  return(-(2/(sigma^2))*V_ring(x))
}
ring_gradlogpdf <- function(x,sigma){
  xgrad = 4*x[1]*(x[1]^2+x[2]^2-1)
  ygrad = 4*x[2]*(x[1]^2+x[2]^2-1)
  output = -(2/(sigma^2))*c(xgrad, ygrad)
  return(output)
}

V_doublewell <- function(x,r){
  if(x>=4){return(6*(x^2)-60)}
  if((x>=0)&(x<4)){return(0.25*(x^4)-2*x^2+4)}
  if((x>=-4/r)&(x<0)){return(0.25*(r*x)^4-2*(r*x)^2+4)}
  if(x<-4/r){return(6*(r*x)^2-60)}
}
grad_V_doublewell <- function(x,r){
  if(x>=4){return(12*(x))}
  if((x>=0)&(x<4)){return(x^3-4*x)}
  if((x>=-4/r)&(x<0)){return(r*(r*x)^3-4*r*(r*x))}
  if(x<-4/r){return(12*r*(r*x))}
}
doublewell_logpdf <- function(x, r, sigma){
  return(-(2/(sigma^2))*V_doublewell(x, r))
}
doublewell_gradlogpdf <- function(x,r, sigma){
  return(-(2/(sigma^2))*grad_V_doublewell(x, r))
}



