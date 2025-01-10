#-----------------------------------
# Likelihood function for BBRM model
#-----------------------------------

msce3_nlm_loglik_IP <- function(par)
{ 
  alpha0 <- 1 # alpha is set to 1 by numberical reasons
  alpha <- matrix(alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  
  
  X0 <- exp(par[1])
  Y1 <- exp(par[2])
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  alpha1 <- par[7]
  gamc <- par[8]
  
  del0 <- exp(-12) # Second mutation rate
  M <-  exp(-13.82) #detection rate 
  alphac <- matrix(alpha1, nrow = gb.nRats, ncol = gb.nTimesteps)
  
  # biological rates
  Nnu0 <- matrix(X0*alpha0,nrow=gb.nRats,ncol=gb.nTimesteps) + X0*alpha0*Y1*df$drate*exp(-Y2*df$drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  # detection rate
  nu2 <- matrix(M, nrow = gb.nRats, ncol = gb.nTimesteps) # M/alpha1
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow=gb.nRats,ncol=gb.nTimesteps) + gam2*(1-exp(-gam1/gam2*df$drate))*expInd
  gammac <- matrix(gamc, nrow = gb.nRats, ncol = gb.nTimesteps)
  
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  
  result <- msce_numerical(t, parList, innerSteps = gb.nInnerSteps)
  haz <- result$hazard
  lnS <- result$lnSurvival
  
  indlik <- (1-df$tumor)*lnS + df$tumor*log(haz*exp(lnS))
  
  return (-sum(indlik))
}

msce3_mle2_loglik_IP <- function(X0,Y1,Y2,gam0,gam1,gam2,alpha1,gamc)
{
  alpha0 <- 1 # alpha is set to 1 by numerical reasons
  alpha <- matrix(alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  
  X0 <- exp(X0)
  alphac <- matrix(alpha1, nrow = gb.nRats, ncol = gb.nTimesteps)
  del0 <- exp(-12)#Second mutation rate
  M <- exp(-13.82) #Detection rate
  Y1 <- exp(Y1)
  
  # biological rates
  Nnu0 <- matrix(X0*alpha0,nrow=gb.nRats,ncol=gb.nTimesteps) +X0*alpha0*Y1*df$drate*exp(-Y2*df$drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  
  # detection rate
  nu2 <- matrix(M, nrow = gb.nRats, ncol = gb.nTimesteps) # M/alpha1
  
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow=gb.nRats,ncol=gb.nTimesteps) + gam2*(1-exp(-gam1/gam2*df$drate))*expInd
  gammac <- matrix(gamc, nrow = gb.nRats, ncol = gb.nTimesteps)
  
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  
  result <- msce_numerical(t,parList, innerSteps = gb.nInnerSteps)
  haz <- result$hazard
  lnS <- result$lnSurvival
  
  indlik <- (1-df$tumor)*lnS + df$tumor*log(haz*exp(lnS))
  
  return (-sum(indlik))
}

msce3_deviance_IP <- function(X0,Y1,Y2,gam0,gam1,gam2,alpha1, gamc, df)
{
  alpha0 <- 1 # alpha is assumed to be 1 by numerical reasons
  alpha <- matrix(alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  
  X0 <- exp(X0)
  Y1 <- exp(Y1)
  alphac <- matrix(alpha1, nrow = gb.nRats, ncol = gb.nTimesteps)
  
  del0 <- exp(-12)#Second mutation rate
  M <- exp(-13.82)#Detection rate
  
  # biological rates
  Nnu0 <- matrix(X0*alpha0,nrow=gb.nRats,ncol=gb.nTimesteps) + X0*alpha0*Y1*df$drate*exp(-Y2*df$drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  
  # detection rate
  nu2 <- matrix(M, nrow = gb.nRats, ncol = gb.nTimesteps)# M/alpha1
  
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow=gb.nRats,ncol=gb.nTimesteps) + gam2*(1-exp(-gam1/gam2*df$drate))*expInd
  gammac <- matrix(gamc, nrow = gb.nRats, ncol = gb.nTimesteps)
  
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  
  result <- msce_numerical(t,parList, innerSteps = gb.nInnerSteps)
  haz <- result$hazard
  lnS <- result$lnSurvival
  
  indlik <- (1-df$tumor)*lnS + df$tumor*log(haz*exp(lnS))
  
  return (-2*sum(indlik))
}

