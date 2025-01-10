#---------------------------------------------------
# Likelihood function for ERR model
#---------------------------------------------------

A_integral <- function(p, tvar){
  result <- 1/(p+2)*(tvar^(p+2) - t[,1]^(p+2)) - t[,1]/(p+1)*(tvar^(p+1) - t[,1]^(p+1))
  return(result)
}

ERR_nlm_loglik <- function(par){
  
  b <- par[1]
  p <- par[2]
  beta1 <- par[3]
  alpha1 <- par[4]
  epsilon <- par[5]
  psi <- par[6]
  tlag <- par[7]
  
  time.matrix <- t
  time.matrix[,3] <- pmax(time.matrix[,3]- tlag, 1)

  acen <- matrix(gb.acen, nrow = gb.nRats, ncol = 1)
  tcen <- matrix(gb.tcen, nrow = gb.nRats, ncol = 1)
  dcen <- matrix(gb.dcen, nrow = gb.nRats, ncol = 1)
  lcen <- matrix(gb.lcen, nrow = gb.nRats, ncol = 1)
  
  k_constant <- beta1*exp(epsilon*(time.matrix[,3] - (t[,1] + (t[,2] - t[,1])/2) - tcen) + psi*(df$drate - dcen))
  
  haz <- exp(b + p*log(t[,3]/lcen))*(1 + k_constant*df$drate*(t[,2] - t[,1]))
  
  lnS <- - 1/p*exp(b +p*log(time.matrix[,3]/lcen)) - exp(b)/p*k_constant*df$drate*(t[,2] - lcen)*exp(p*(time.matrix[,3]-lcen)) + exp(b)/(p^2)*k_constant*df$drate*exp(p*(time.matrix[,3]-lcen)) - exp(b)/(p^2)*k_constant*df$drate*exp(p*(t[,1]-lcen))
  
  indlik <- (1-df$tumor)*lnS + df$tumor*(log(haz) + lnS)
  return (-sum(indlik))
}

ERR_mle_loglik <- function(f, p, beta1, alpha1, epsilon, psi, tlag){
  
  acen <- matrix(gb.acen, nrow = gb.nRats, ncol = 1)
  tcen <- matrix(gb.tcen, nrow = gb.nRats, ncol = 1)
  dcen <- matrix(gb.dcen, nrow = gb.nRats, ncol = 1)
  lcen <- matrix(gb.lcen, nrow = gb.nRats, ncol = 1)
  
  time.matrix <- t
  time.matrix[,3] <- pmax(time.matrix[,3]- tlag, 1)


  k_constant <- beta1*exp(epsilon*(time.matrix[,3] - (t[,1] + (t[,2] - t[,1])/2) - tcen) + psi*(df$drate - dcen))
  
  haz <- exp(b + p*log(time.matrix[,3]/lcen))*(1 + k_constant*df$drate*(t[,2] - t[,1]))
  
  lnS <- - 1/p*exp(b +p*log(time.matrix[,3]/lcen)) - exp(b)/p*k_constant*df$drate*(t[,2] - lcen)*exp(p*(t[,3]-lcen)) + exp(b)/(p^2)*k_constant*df$drate*exp(p*(time.matrix[,3]-lcen)) - exp(b)/(p^2)*k_constant*df$drate*exp(p*(t[,1]-lcen))
  
  indlik <- (1-df$tumor)*lnS + df$tumor*(log(haz) + lnS)

  return (-sum(indlik))
}

ERR_deviance<- function(f, p, beta1, alpha1, epsilon, psi){
  acen <- matrix(gb.acen, nrow = gb.nRats, ncol = 1)
  tcen <- matrix(gb.tcen, nrow = gb.nRats, ncol = 1)
  dcen <- matrix(gb.dcen, nrow = gb.nRats, ncol = 1)
  lcen <- matrix(gb.lcen, nrow = gb.nRats, ncol = 1)
  tlag <- par[7]
  
  time.matrix <- t
  time.matrix[,3] <- pmax(time.matrix[,3]- tlag, 1)

  k_constant <- beta1*exp(epsilon*(time.matrix[,3] - (t[,1] + (t[,2] - t[,1])/2) - tcen) + psi*(df$drate - dcen))
  
  haz <- exp(b + p*log(time.matrix[,3]/lcen))*(1 + k_constant*df$drate*(t[,2] - t[,1]))
  
  lnS <- - 1/p*exp(b +p*log(time.matrix[,3]/lcen)) - exp(b)/p*k_constant*df$drate*(t[,2] - lcen)*exp(p*(t[,3]-lcen)) + exp(b)/(p^2)*k_constant*df$drate*exp(p*(time.matrix[,3]-lcen)) - exp(b)/(p^2)*k_constant*df$drate*exp(p*(t[,1]-lcen))
  
  indlik <- (1-df$tumor)*lnS + df$tumor*(log(haz) + lnS)

  return (-2*sum(indlik))
  
  
}
