#-----------------------------
# Code for computing relevant times and time line of carcinogenesis
# used in paper
# ST September 2024
#----------------------------
library(msce)
library(ggplot2)
# All time information is considered in weeks and in matrix form compatible
# with msce functions. Biological parameters are given by 
# (logY0, logY1, Y2 = 0, gamma_0, gamma_1, gamma_2, logG, Gamma)
wparall <- c(2.14, exp(5.02e-1), 0, 3.98e-2, 416e-3, 0.416e-1, exp(1.98), 1.39e-1)
wparadeno <- c(6.6e-1, exp(5.15e-1), 0, 4.94e-2, 1.07e-3, 0.94e-1, exp(1.6), 1.57e-1 )
wparsq <- c(2.8e-1, 0, 0, 4.97e-2, 2.13e-3, 3.27e-1, exp(1.6), 2.21e-1)
#------------------------------
# Timeline of carcinogenesis
#------------------------------

# Initiation period

# Survival for no TIC and no CSC

Surv_0 <- function(t, par, drate, tb, expdur) {
  # par are the biological parameters and should be given in the form above
  # the exposure rate drate is in WL
  # tb is the age at begin of exposure in weeks and expdur is exposure duration
  
  # compute the relevant biological parameters
  Nnu0 <- exp(par[1])
  Nnu1 <- exp(par[1]) + exp(par[1])*exp(par[2])*drate*exp(-par[3]*drate)
  Nnu2 <- Nnu0 # first mutation rate return to level pre-exposition
  
  # booleans for which interval and corresponding function
  before<- t < tb
  after <- t > tb + expdur
  during <- t < tb + expdur & t > tb
  exponent_0 <- -Nnu0*t*before + during*(-Nnu1*t - (Nnu0 - Nnu1)*tb) + after*(-Nnu0*t + expdur*(Nnu0 - Nnu1))
  S_0 <- exp(exponent_0)
  
}

# Compute initiation period: integral of Surv_0 = \Delta_I = T_I
T_I <- function(par, drate, tb, expdur) {
  Surv_to_integrate <- function(t) {
    return(Surv_0(t, par, drate, tb, expdur))
  }
  integrate(Surv_to_integrate, 2/7, 300)$value # for numerical reasons take finite integral 
}

#------------------------------------------------------------
# Hyperplasia formation period

# Survival for no CSC
Surv_C <- function(t, par, drate, tb, expdur){
  
  # time information in weeks!
  
  Y0 <- exp(par[1])
  Y1 <- par[2]
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  
  # constant biological parameter
  alpha0 <- matrix(1, nrow = length(t), ncol = 3) #TIC division rate
  del0 <- exp(-12)#10^-6
  
  # produces the time matrix in msce format
  time_matrix <- matrix(0, nrow = length(t), ncol = 3)
  before <- t <= tb
  during <- tb < t & t <= tb + expdur
  after <- t > tb + expdur 
  time_matrix[, 1] <- after*tb
  time_matrix[, 2] <- after*(tb + expdur) + during*tb
  time_matrix[, 3] <- t
  
  expInd <- matrix(0, nrow = length(t), ncol = 3)
  expInd[ , 2] <- after
  expInd[ , 3] <- during
  
  # biological rates
  Nnu0 <- matrix(Y0*alpha0,nrow=length(t),ncol=3) + Y0*alpha0*Y1*drate*exp(-Y2*drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=length(t),ncol=3)
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow= length(t),ncol=3) + gam2*(1-exp(-gam1/gam2*drate))*expInd
  
  #use TSCE
  parList <- list(Nnu0=Nnu0, alpha=alpha0, gamma=gamma, nu1=nu1)
  result_tsce <- tsce(time_matrix,parList)
  lnS_tsce <- result_tsce$lnSurvival
  S_tsce <- exp(lnS_tsce)
  
  return(S_tsce)
  
}

# T_C: mean age at first CSC is the integral of the survival function Surv_C
T_C <- function(par, drate, tb, expdur) {
  Surv_to_integrate <- function(t) {
    return(Surv_C(t, par, drate, tb, expdur))
  }
  integrate(Surv_to_integrate, 2/7, 300)$value # for numerical reasons take finite integral 
}

# Difference of survival functions S_C - S_0 
Delta_S_C0 <- function(t, par, drate, tb, expdur){
  return(Surv_C(t, par, drate, tb, expdur) - Surv_0(t, par, drate, tb, expdur))
}

# Hyperplasia formation period
Delta_H <- function( par, drate, tb, expdur){
  Surv_to_integrate <- function(t) {
    return(Delta_S_C0(t, par, drate, tb, expdur))
  }
  integrate(Surv_to_integrate, 2/7, 300)$value # for numerical reasons take finite integral 
}
  
  
# Normalized difference needed for Figure 7
Delta_S_C0_normalized <- function(t, par, drate, tb, expdur) {
  Delta_denominator <- Delta_H(par, drate, tb, expdur)
  unnorm_surv <- Delta_S_C0(t, par, drate, tb, expdur)
  return(unnorm_surv/Delta_denominator)
}

#---------------------------------------------------
# Tumor formation period

# Survival for no detected tumor

Surv_D <- function(t, par, drate, tb, expdur){
  
  X0 <- exp(par[1])
  Y1 <- par[2]
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  alpha1 <- par[7]
  gamc <- par[8]
  
  # constant biological parameters
  alpha0 <- matrix(1, nrow = length(t), ncol = 3) #gb.alpha0
  del0 <- exp(-12)
  M <-  exp(-13.82) 
  alphac <- matrix(alpha1, nrow = length(t), ncol = 3)
  
  # produces the time matrix in msce format
  # produces the time matrix in msce format
  time_matrix <- matrix(0, nrow = length(t), ncol = 3)
  before <- t <= tb
  during <- tb < t & t <= tb + expdur
  after <- t > tb + expdur 
  time_matrix[, 1] <- after*tb
  time_matrix[, 2] <- after*(tb + expdur) + during*tb
  time_matrix[, 3] <- t
  
  expInd <- matrix(0, nrow = length(t), ncol = 3)
  expInd[ , 2] <- after
  expInd[ , 3] <- during
  
  # biological rates
  Nnu0 <- matrix(X0*alpha0,nrow=length(t),ncol=3) + X0*alpha0*Y1*drate*exp(-Y2*drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=length(t),ncol=3)
  # detection rate
  nu2 <- matrix(M, nrow = length(t), ncol = 3) # M/alpha1
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow= length(t),ncol=3) + gam2*(1-exp(-gam1/gam2*drate))*expInd
  gammac <- matrix(gamc, nrow = length(t), ncol = 3)
  
  # use MSCE
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha0,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  result_msce <- msce_numerical(time_matrix, parList, innerSteps = 1000)
  lnS_msce <- result_msce$lnSurvival
  S_msce <- exp(lnS_msce)
  
  return(S_msce)
}

# Mean age at diagnosis is the integral of survival for no detected tumor
T_D <- function(par, drate, tb, expdur) {
  Surv_to_integrate <- function(t) {
    return(Surv_D(t, par, drate, tb, expdur))
  }
  integrate(Surv_to_integrate, 2/7, 300)$value # for numerical reasons take finite integral 
}

# Difference of survival functions Surv_D - Surv_C 
Delta_S_DC <- function(t, par, drate, tb, expdur){
  return(Surv_D(t, par, drate, tb, expdur) - Surv_C(t, par, drate, tb, expdur))
}

# Tumor formation period
Delta_T <- function( par, drate, tb, expdur){
  Surv_to_integrate <- function(t) {
    return(Delta_S_DC(t, par, drate, tb, expdur))
  }
  integrate(Surv_to_integrate, 2/7, 300)$value # for numerical reasons take finite integral 
}


# Normalized difference needed for Figure 7
Delta_S_DC_normalized <- function(t, par, drate, tb, expdur) {
  Delta_denominator <- Delta_T(par, drate, tb, expdur)
  unnorm_surv <- Delta_S_DC(t, par, drate, tb, expdur)
  return(unnorm_surv/Delta_denominator)
}


#----------------------------------
# Hazard of detected cancer
#----------------------------------
Haz_D <- function(t, par, drate, tb, expdur){
  
  X0 <- exp(par[1])
  Y1 <- par[2]
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  alpha1 <- par[7]
  gamc <- par[8]
  
  # constant biological parameters
  alpha0 <- matrix(1, nrow = length(t), ncol = 3) #gb.alpha0
  del0 <- exp(-12)
  M <-  exp(-13.82) 
  alphac <- matrix(alpha1, nrow = length(t), ncol = 3)
  
  # produces the time matrix in msce format
  # produces the time matrix in msce format
  time_matrix <- matrix(0, nrow = length(t), ncol = 3)
  before <- t <= tb
  during <- tb < t & t <= tb + expdur
  after <- t > tb + expdur 
  time_matrix[, 1] <- after*tb
  time_matrix[, 2] <- after*(tb + expdur) + during*tb
  time_matrix[, 3] <- t
  
  expInd <- matrix(0, nrow = length(t), ncol = 3)
  expInd[ , 2] <- after
  expInd[ , 3] <- during
  
  # biological rates
  Nnu0 <- matrix(X0*alpha0,nrow=length(t),ncol=3) + X0*alpha0*Y1*drate*exp(-Y2*drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=length(t),ncol=3)
  # detection rate
  nu2 <- matrix(M, nrow = length(t), ncol = 3) # M/alpha1
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow= length(t),ncol=3) + gam2*(1-exp(-gam1/gam2*drate))*expInd
  gammac <- matrix(gamc, nrow = length(t), ncol = 3)
  
  # use MSCE
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha0,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  result_msce <- msce_numerical(time_matrix, parList, innerSteps = 1000)
  haz_msce <- result_msce$haz
  
  return(haz_msce)
}


#--------------------------------------------
# normalized probability of cancer detection

# Normalized probability of tumor detection
Surv_Haz_D <- function(t, par, drate, tb, expdur) {
  norm_surv <- Surv_D(t, par, drate, tb, expdur)*Haz_D(t, par, drate, tb, expdur)
  # already normalized since haz = S'/S
  return(norm_surv)
}


#-------------------
# Figure 7
#-------------------

# The following code produces the bottom right panel of figure 7
cbPalette <- c("#999999",  "#009E73", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot(data.frame(x = c(2/7,250)), aes(x = x))+ 
  geom_path(aes(colour= cbPalette[3]), stat="function", fun = Delta_S_C0_normalized, args = list(par =wparadeno, drate = 60, tb = 13, expdur = 51), size = 1.1)+
  geom_path(aes(colour= cbPalette[5]), stat="function", fun = Delta_S_DC_normalized , args = list(par =wparadeno, drate = 60, tb = 13, expdur = 51), size = 1.1)+
  geom_path(aes(colour= cbPalette[6]), stat="function", fun = Surv_Haz_D, args = list(par =wparadeno, drate = 70, tb = 13, expdur = 0.6), size = 1.1)+
  scale_colour_identity(" ", guide="legend", 
                        labels = c( "hyperplasia formation", "tumor formation", "tumor detection"), 
                        breaks = cbPalette[c(3,5,6)]) +
  xlab("Rat age (weeks)") +
  ylab("Probability distribution") +
  #theme(legend.text = element_text(
    #face = "bold",
   # size = 14)
  #) +
  theme(legend.position = c(.75,.85)) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(axis.title.y = element_text( size = 14))
p

#---------------------------------
# Sojourn time A_H
#---------------------------------

# unconditioned clone survival
surv_clone <- function(t, par, drate, t_b, expdur){
  t_b <- t_b
  gam0 <- par[4] #
  gam1 <- par[4] + par[6]*(1- exp(-par[5]/par[6]*drate)) #
  gam2 <- gam0
  

  Gs <- c(1,1,1) # we assumed alpha constant = 1
  Ds <- c(1-gam0, 1-gam1, 1-gam2)
  Ms <- c(exp(-12), exp(-12), exp(-12)) # second mutation rates
  Deltas <- sqrt((Ds+Gs+Ms)^2 - 4*Ds*Gs)
  
  As <- (Ms + Ds + Gs+ Deltas)/(2*Gs)
  Bs <- (Ms + Ds + Gs- Deltas)/(2*Gs)
  
  A0 <- As[1]
  Ae <- As[2]
  Aa <- As[3]
  
  B0 <- Bs[1]
  Be <- Bs[2]
  Ba <- Bs[3]
  
  # p_inf to make the probability go to 1
  p_num <- (1-B0)*exp((A0-B0)*t_b)*((Be-B0)*(A0-Ae) + (A0 -Be)*(B0 -Ae)*exp((Ae-Be)*expdur)) + (1-A0)*(Be-B0)*(Ae-B0)*(1 - exp(Ae-Be)*expdur)
  p_den <- exp((A0-B0)*t_b)*((Be-B0)*(A0-Ae) + (A0 -Be)*(B0 -Ae)*exp((Ae-Be)*expdur)) + (Be-B0)*(Ae-B0)*(1 - exp((Ae-Be)*expdur))
  p_inf <- p_num/p_den
  
  
  # Survival before exposure
  fun1 <- function(s) {
    A <- As[1]
    B <- Bs[1]
    c <- (A -1)/(1 - B)
    St <- (A + B*c*exp((A-B)*s))/(1 + c*exp((A-B)*s)  )
    prob <- 1-St
    prob <- prob/p_inf
    return(1 - prob)
  }

  # Survival during exposure
  fun2 <- function(s) {
    A0 <- As[1]
    B0<- Bs[1]
    A <- As[2]
    B <- Bs[2]
    Numerator <- exp((A-B)*(s-t_b))*(1 -A)*(A*(B-B0) + B0*exp((A0-B0)*t_b)*(A0 - B)) + (B -1)*(A0*(A-B0) + B0*exp((A0-B0)*t_b)*(A0 -A) )
    Denominator <- exp((A-B)*(s-t_b))*(1 -A)*((B-B0) + exp((A0-B0)*t_b)*(A0 - B)) + (B -1)*((A-B0) + exp((A0-B0)*t_b)*(A0 -A) )
    
    St <- Numerator/Denominator
    
    prob <- 1-St
    prob <- prob/p_inf
    return(1 - prob)
  }
  
  # After exposure
  fun3 <- function(s) {
    
    Numerator <- (1-Aa)*exp((Aa-Ba)*(s-t_b-expdur ))*(B0*exp((A0-B0)*t_b)*((Be-B0)*(A0-Ae) + (A0-Be)*(B0-Ae)*exp((Ae-Be)*expdur)) + A0*(Be-B0)*(Ae-B0)*(1 - exp((Ae-Be)*expdur)))
    Numerator <- Numerator + (B0 -1)*(B0*exp((A0-B0)*t_b)*(A0-Ae)*(Be-A0)*(1-exp((Ae-Be)*expdur)) + A0*((A0-Ae)*(Be-B0)*exp((Ae-Be)*expdur) + (Be-A0)*(Ae-B0)) )
    
    Denominator <- (1-Aa)*exp((Aa-Ba)*(s-t_b-expdur ))*(exp((A0-B0)*t_b)*((Be-B0)*(A0-Ae) + (A0-Be)*(B0-Ae)*exp((Ae-Be)*expdur)) + (Be-B0)*(Ae-B0)*(1 - exp((Ae-Be)*expdur)))
    Denominator <- Denominator + (B0 -1)*(exp((A0-B0)*t_b)*(A0-Ae)*(Be-A0)*(1-exp((Ae-Be)*expdur)) + ((A0-Ae)*(Be-B0)*exp((Ae-Be)*expdur) + (Be-A0)*(Ae-B0)) )
    
    St <- Numerator/Denominator
    
    prob <- 1-St
    prob <- prob/p_inf
    return(1 - prob)
  }
  # boolean to locate t in the intervals
  before <- t <= t_b
  during <- t_b < t & t <= t_b + expdur
  after <- t_b + expdur < t
  
  retval <- before*fun1(t) + during*fun2(t) + after*fun3(t)
  
  return(retval)
  
}

compute_soj <- function(par, drate, expdur){
  surv_to_integrate <- function(s) {
    return(surv_clone(t = s, par = par, drate = drate, t_b = 13, expdur = expdur))
  }
  integrate(surv_to_integrate, 0, 300)$value
  
}

compute_soj(wparall, 0, 0)



# clone survival conditioned on rat survival

surv_clone_cond_function <- function(t, par, drate, t_b, expdur){

  to_integrate <- function(s) {
    num <- surv_clone(t, par, drate, t_b, expdur) - surv_clone(s, par, drate, t_b, expdur)
    denom <- 1 - surv_clone(s, par, drate, t_b, expdur) 
    surv_rat <- Surv_D(s, par, drate, tb =  t_b, expdur)
    haz_rat <- Haz_D(s, par, drate, tb = t_b, expdur)
		
    return(surv_rat*haz_rat*num/denom)
  }

  retval <- sapply(t, function (t) {integrate(to_integrate, t, 300)$value  })

  return(retval)
}


#conditioned sojourn time

compute_soj_cond <- function(par, drate, expdur){

  surv_clone_cond <- function(t, par, drate, t_b, expdur){
    to_integrate <- function(s) {
      surv_clone <- surv_clone(s, par, drate, t_b, expdur)		
      return(surv_clone)
    }

  S_c <- surv_clone(t, par, drate, t_b, expdur)
  denom <- 1 - S_c
  num <- sapply(t, function (t) {integrate(to_integrate, 0, t)$value  })
  num <- num - S_c*t

  retval <- num/denom
  return(retval)
  }

  surv_to_integrate <- function(s) {
    clone <-surv_clone_cond(s, par, drate, t_b = 13, expdur)
    surv_rat <- Surv_D(s, par, drate, tb =  13, expdur)
    haz_rat <- Haz_D(s, par, drate, tb = 13, expdur)
    return(clone*surv_rat*haz_rat)

  }
  integrate(surv_to_integrate, 1, 300)$value
  
}
compute_soj_cond(wparall, 0, 0)

#--------------------------------------------------------------
# Produce Table 4: conditioned and unconditioned sojourn times
#--------------------------------------------------------------


# Response to different exposure scenarios
soj_exp <- data.frame(row.names = c("All_0", "All_Ll", "All_Ml","All_Mh", "All_Hl", "All_Hh", "Adeno_0" , "Adeno_Ll", "Adeno_Ml", "Adeno_Mh", "Adeno_Hl", "Adeno_Hh", "Squamous_0" , "Squamous_Ll", "Squamous_Ml", "Squamous_Mh", "Squamous_Hl", "Squamous_Hh"), Soj_uncond = numeric(18), Soj_cond = numeric(18), T_D = numeric(18))

exp_scen <- data.frame(row.names = c("0", "Ll", "Ml", "Mh", "Hl", "Hh"), Dose_rate = c(0, 70, 52, 565, 60, 390), Dose = c(0, 42, 416, 565, 3060, 5460))

list_par <- list(All = wparall, Adeno = wparadeno, Squamous = wparsq)

soj_exp

for (hist in 1:3){
 	upar <- list_par[[hist]]
  
	for (iter in 1:6){
		drate <- exp_scen$Dose_rate[iter]
		expdur <- exp_scen$Dose[iter]/drate

		if (drate == 0){
			expdur <- 0
		}
    
		soj_exp$Soj_uncond[6*(hist-1) + iter] <- compute_soj(upar, drate, expdur)
		soj_exp$Soj_cond[6*(hist-1) + iter] <- compute_soj_cond(upar, drate, expdur)
		soj_exp$T_D[6*(hist-1) + iter]<- T_D(upar, drate, tb = 13, expdur)
    
	}
}


soj_exp


soj_exp <- round(soj_exp)
soj_exp

soj_exp$Check <- soj_exp$T_D > soj_exp$Soj_cond
soj_exp

write.csv(soj_exp, file = "I:/Rats/Drafts/Soj_cond.csv")












#-------------------------------
# Produce Fig S6
ages <- c(0:250)

SMh <-  numeric(251)
SHl <- numeric(251)
SHh <- numeric(251)

SMh[1] <- 1
SHl[1] <- 1
SHh[1] <- 1

for (iter in 1: 250) {
  SMh[iter+1] <- surv_clone_cond_function(ages[iter +1], wparadeno, drate = 565, t_b = 13, expdur = 1)
  SHl[iter+1] <- surv_clone_cond_function(ages[iter +1], wparadeno, drate = 60, t_b = 13, expdur = 51)
  SHh[iter+1] <- surv_clone_cond_function(ages[iter +1], wparadeno, drate = 390, t_b = 13, expdur = 14)
}

hDd <- data.frame(
  x = rep(ages, 1),
  y = c(SMh, SHl, SHh),
  expgrp = rep( c("Mh","Hl","Hh" ), each = length(ages)),
  cond = "Conditioned"
)
head(hDd)


SMhu <-  numeric(251)
SHlu <- numeric(251)
SHhu <- numeric(251)

SMhu[1] <- 1
SHlu[1] <- 1
SHhu[1] <- 1


for (iter in 1: 250) {
  SMhu[iter+1] <- surv_clone(ages[iter+1], wparadeno, drate = 565, t_b = 13, expdur = 1)
  SHlu[iter+1] <- surv_clone(ages[iter+1], wparadeno, drate = 60, t_b = 13, expdur = 51)
  SHhu[iter+1] <- surv_clone(ages[iter+1], wparadeno, drate = 390, t_b = 13, expdur = 14)
}

hDdu <- data.frame(
  x = rep(ages, 1),
  y = c(SMhu, SHlu, SHhu),
  expgrp = rep( c("Mh","Hl","Hh" ), each = length(ages)),
  cond = "Unconditioned"
)
pf.md <- rbind(hDd,hDdu)
names(pf.md) <- c("x", "y", "expgrp", "Model")
head(pf.md)


ggplot( data = pf.md) + 
  geom_line(aes(x=x, y=y, colour = expgrp, linetype = Model), size =1) +
  coord_cartesian(xlim = c(0, 175), expand = FALSE)+
  scale_fill_manual(name = "Exposure class", breaks =c("Mh","Hl","Hh" ), values= cbPalette[4:6], aesthetics = c("color", "fill"))+
  scale_linetype_manual(name = " ", breaks = c("Conditioned", "Unconditioned"), values = c(1,2)) +
  xlab("Rat age (weeks)") + ylab("Surviving fraction of Hyperplasia")+
  theme(text = element_text(size=12),legend.position = c(.9,0.75)) 


ggsave(file="ADC_conduncond.eps")


#------------------------------------------------------------------------
# Produce Fig S6: dependence of clone survival on time of origin of clone

wparadeno <- c(6.6e-1, 5.15e-1, 0, 4.94e-2, 1.07e-3, 0.94e-1, exp(1.6), 1.57e-1 )

# The following code produces the bottom right panel of figure 7
cbPalette <- c("#999999",  "#009E73", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot(data.frame(x = c(2/7,250)), aes(x = x))+ 
  geom_path(aes(colour= cbPalette[1]), stat="function", fun = compute_soj_time, args = list(par =wparadeno, drate = 0, expdur = 0), size = 1.1)+
  geom_path(aes(colour= cbPalette[2]), stat="function", fun = compute_soj_time, args = list(par =wparadeno, drate = 70, expdur = 0.6), size = 1.1)+
  geom_path(aes(colour= cbPalette[3]), stat="function", fun = compute_soj_time, args = list(par =wparadeno, drate = 52, expdur = 8), size = 1.1)+
  geom_path(aes(colour= cbPalette[4]), stat="function", fun = compute_soj_time, args = list(par =wparadeno, drate = 565, expdur = 1), size = 1.1)+
  geom_path(aes(colour= cbPalette[5]), stat="function", fun = compute_soj_time, args = list(par =wparadeno, drate = 60, expdur = 51), size = 1.1)+
  geom_path(aes(colour= cbPalette[6]), stat="function", fun = compute_soj_time, args = list(par =wparadeno, drate = 390, expdur = 14), size = 1.1)+

  scale_colour_identity(" ", guide="legend", 
                        labels = c( "Unexposed", "Ll", "Ml", "Mh", "Hl", "Hh"), 
                        breaks = cbPalette) +
  xlab("") +
  ylab("Probability distribution") +
  #theme(legend.text = element_text(
    #face = "bold",
   # size = 14)
  #) +
  theme(legend.position = c(.75,.85)) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(axis.title.y = element_text( size = 14))
p


#----------------------------------------------
# Produce Fig S8

# Mean sojourn time as function of age at end of follow up (Eq S.26)

compute_soj_fixed_td <- function(t_d, par, drate, t_b, expdur){
  surv_to_integrate <- function(s) {
    numerator <- surv_clone(s, par, drate, t_b,  expdur ) - surv_clone(t_d, par, drate, t_b,  expdur )
    denominator <- 1 - surv_clone(t_d, par, drate, t_b,  expdur )
    return(numerator/denominator)
  }
  integrate(surv_to_integrate, 0, t_D)
  
}

# The following code produces the bottom right panel of figure 7
cbPalette <- c("#999999",  "#009E73", "#56B4E9", "#0072B2", "#D55E00", "#CC79A7")
p <- ggplot(data.frame(x = c(2/7,250)), aes(x = x))+ 
  geom_path(aes(colour= cbPalette[1]), stat="function", fun = compute_soj_fixed_td, args = list(par =wparadeno, drate = 0, t_b = 13, expdur = 0), size = 1.1)+
  geom_path(aes(colour= cbPalette[2]), stat="function", fun = compute_soj_fixed_td, args = list(par =wparadeno, drate = 70, t_b = 13, expdur = 0.6), size = 1.1)+
  geom_path(aes(colour= cbPalette[3]), stat="function", fun = compute_soj_fixed_td, args = list(par =wparadeno, drate = 52, t_b = 13, expdur = 8), size = 1.1)+
  geom_path(aes(colour= cbPalette[4]), stat="function", fun = compute_soj_fixed_td, args = list(par =wparadeno, drate = 565, t_b = 13, expdur = 1), size = 1.1)+
  geom_path(aes(colour= cbPalette[5]), stat="function", fun = compute_soj_fixed_td, args = list(par =wparadeno, drate = 60, t_b = 13, expdur = 51), size = 1.1)+
  geom_path(aes(colour= cbPalette[6]), stat="function", fun = compute_soj_fixed_td, args = list(par =wparadeno, drate = 390, t_b = 13, expdur = 14), size = 1.1)+

  scale_colour_identity(" ", guide="legend", 
                        labels = c( "Unexposed", "Ll", "Ml", "Mh", "Hl", "Hh"), 
                        breaks = cbPalette) +
  xlab("") +
  ylab("Probability distribution") +
  #theme(legend.text = element_text(
    #face = "bold",
   # size = 14)
  #) +
  theme(legend.position = c(.75,.85)) +
  theme(axis.title.x = element_text( size = 14)) +
  theme(axis.title.y = element_text( size = 14))
p




#-----------------------------------
# Mean tumor age A_T
#-----------------------------------
mean_tumAge <- function(alpha1,gamc){
  M <- exp(-13.82)# constant detection rate
  G <- alpha1 # CSCs duplication rate
  A_T <- abs(log(gamc^2/(G*M))/gamc)
  return (A_T)
}
