#--------------------------------------------------------------------------
# Fit parameters of ERR model for PNNL rats
#--------------------------------------------------------------------------

rm(list = ls()) # remove all objects from the current workspace

library(forcats)
library(ggplot2)
library(scales)
library(gplots)
library(survival)
library(survminer)
library(cowplot)

library(msce)
library(gnm) # generalized non-linear models
library(bbmle) # mle2 fitting with functions similar to glm fitting
#library(splines)
library(numDeriv) # for hessian(.)
library(matrixcalc) # for testing if hessian is positive definite
library(MASS)
library(Formula)

# http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/
# The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#------------------------------------------------------------------
# read prepared data filtered from rat.wh
# from Heidenreich et al. RasRes 151, 209-17, 1999
#------------------------------------------------------------------
setwd("I:/Rats/stats/data")
dataSet <- c("PNNLrats-231121.Rdata")
load(file = dataSet)
setwd("I:/Rats/stats")

dim(df.pnnl) # 4260   24
names(df.pnnl)
# [1] "id"        "aboexp"    "aeoexp"    "aeofup"    "expdur"    "agrp"      "exposed"   "dose"      "drate"     "hWLM"     
#[11] "hWL"       "clsWLM"    "clsWLM4"   "clsWL"     "clsWL3"    "tumor"     "tumfat"    "adeno"     "epid"      "adesq"    
#[21] "sarco"     "deathcode" "dae"       "dde"  

#-------------------------------
# set data frane for analysis
#-------------------------------

df <- df.pnnl
nRats <- dim(df)[1]
nTimesteps <- 3

#-----------------------------------------------------------------------------
# define time matrix in msce format
#-----------------------------------------------------------------------------
# time scaling constants
dpw <- 7 # days per week 
dpm <- 365/12 # days per month

#  matrix with time steps
# define t and set to zero
t <- matrix(0,nrow=nRats,ncol=nTimesteps)

# death after exposure
t[df$dae,1] <- df$aboexp[df$dae]/dpw
t[df$dae,2] <- df$aeoexp[df$dae]/dpw

# death during exposure
t[df$dde,2] <-df$aeoexp[df$dde]/dpw

# end of follow up always comes last
t[,3] <- df$aeofup/dpw

# define exposure index matrix and set to zero
expInd <- matrix(0,nrow=nRats,ncol=nTimesteps)

# death after exposure
expInd[df$dae,2] <- 1

# death during exposure
expInd[df$dde,3] <- 1

#-----------------------
# global variables

gb.acen <- median(t[,1] + (t[,2] - t[,1])/2)
gb.tcen <- median(t[,3] - (t[,1] + (t[,2] - t[,1])/2))
gb.dcen <- median(df$drate)
gb.lcen <- median(t[,3])


# msce_numerical
gb.nInnerSteps <- 1000 # default

nInnerSteps <- 1000 # inner steps for numerical fitting

#-----------------------------------
# select models and likelihoods
#-----------------------------------
df$tum <- df$tumor # boolean for presence of tumor

# model solver
mdsolv <- character()
mdsolv <- "tsce-TG" # msce solver

# uncomment for histological subtypes
histo <- character()
histo <- "all"
#histo <- "adeno"
#histo <- "epid"

# variable tumor used in likelihood function: denotes the histology
 
if (histo == "all") {
	df$tumor <- df$tum
}

if (histo == "adeno"){
	df$tumor <- df$adeno
}
if (histo == "epid") {
	df$tumor <- df$epid
}

#---------------------------------------------------------------------------------
# set data frame for analysis
#---------------------------------------------------------------------------------
df <- df.pnnl
#df <- subset(df.pnnl, exposed == "no")
nRats <- dim(df)[1]
nTimesteps <- 3
source("I:/Rats/stats/subscripts/globvar.R")
#-----------------------------------------------------------------------------
# define matrices in msce format
#-----------------------------------------------------------------------------
# time scaling constants
dpw <- 7 # days per week 
dpm <- 365/12 # days per month

#  matrix with time steps different from MSCE
# define t and set to zero
t <- matrix(0,nrow=nRats,ncol=nTimesteps)

# death after exposure
t[,1] <- df$aboexp/dpw
t[df$dae,2] <- df$aeoexp[df$dae]/dpw

# death during exposure, end of exposure is also end of follow up
t[df$dde,1] <- df$aeoexp[df$dde]/dpw
t[df$dde,2] <-df$aeofup[df$dde]/dpw

# end of follow up always comes last
t[,3] <- df$aeofup/dpw

# no need for expInd bc t is already designed for right cumulative dose

oldInd <- matrix(0, nrow = gb.nRats, ncol = 1)
for (iter in 1: nRats){
  oldInd[iter] <- t[iter, 3 ] > median(t[,3])
}

head(1-oldInd)

# msce_numerical
#nInnersteps <- 1000 # default


gb.acen <- median(t[,1] + (t[,2] - t[,1])/2)
gb.tcen <- median(t[,3] - (t[,1] + (t[,2] - t[,1])/2))
gb.dcen <- median(df$drate)
gb.lcen <- median(t[,3])

#--------------------------------
# import likelihood functions
#--------------------------------
source("I:/Rats/stats/subscripts/ERR-likeli.R")

# controls CI printing
has.profile <- FALSE
has.Hessian <- FALSE

smloglik <- ERR_mle_loglik
snloglik <- ERR_nlm_loglik
sdeviance <- ERR_deviance

#----------------------------
# read estimated parameters
#----------------------------

parDF <- read.csv("I:/Rats/stats/results/ERR/ERR-step/ERR-step-epid2_parms.csv")


setwd("I:/Rats/stats")
upar <- vector()
upar <- as.numeric(parDF$parval[c(1:7)])
upar

wpar <- upar

wpar<- numeric(7)
wpar<- upar

names(wpar) <- c("b", "p", "beta1", "alpha1", "epsilon", "psi", "tlag")
wpar

sdeviance(wpar[1], wpar[2], wpar[3], wpar[4], wpar[5], wpar[6])

#------------------------------------------------------------------
# Fit parameters with bbmle package

mmcall = 0
durat<- system.time(
  mle.3IP <- mle2(minuslog = smloglik,
                  start=list(f = wpar[1], p = wpar[2], beta1= wpar[3], alpha1 = wpar[4], epsilon = wpar[5], psi = wpar[6], tlag = wpar[7]), 
                  parameters=list(f~1, p~1, p_old~1, beta1~1, alpha1~1, epsilon~1, psi~1, tlag~1), 
                  fixed=list(
                    b = wpar[1]
                    #p = wpar[2],
                    #beta1 = wpar[3],
                    #alpha1 = wpar[4],
                    #epsilon = wpar[5],
                    #psi =wpar[6],
			  tlag = wpar[7]
                  ),
                  lower = c( p = 1,  beta1 = 5e-3, epsilon = -1e-1, alpha1 = -7e-2, psi = -5e-3), #  
                  upper = c( p = 10, beta1 = 7e-2, epsilon = -1e-2, alpha1 = - 5e-3, psi = -1e-3 ), #
                  
                  
                  method = "L-BFGS-B",
                  #method = "CG",
                  #method = "Brent",
                  data = df)
) # system.time
durat

cat(sprintf("Total calls: %d", mmcall))

summary(mle.3IP)

wpar <- coef(mle.3IP)
names(wpar) <- c("b", "p", "beta1", "alpha1", "epsilon", "psi")


#------------------------------------
# Survival function for ERR model
#------------------------------------
# constants
df$aeoexp.true <- df$aeoexp
df$aeoexp.true[df$dde] <- df$aeofup[df$dde] # true end of exposure = end of follow up
df$expdur.true <- df$aeoexp.true- df$aboexp # true exposure duration
df$tsme <- df$aeofup - (df$aboexp + (df$aeoexp.true- df$aboexp)/2) # time since median exposure
df$aame <- df$aboexp + (df$aeoexp.true- df$aboexp)/2 # time since median exposure

df$tsme <- df$aeofup - (df$aboexp + (df$aeoexp.true- df$aboexp)/2) # time since median exposure
df$aame <- df$aboexp + (df$aeoexp.true- df$aboexp)/2 # time since median exposure

gb.tcen <- median(df$tsme)/7
gb.dcen <- median(df$drate)
gb.lcen <- median(df$aeofup)/7
gb.acen <- median(df$aeofup)/7
gb.tlag <- 15 # weeks


A_integral <- function(p, tvar){
  result <- 1/(p+2)*(tvar^(p+2) - t[,1]^(p+2)) - t[,1]/(p+1)*(tvar^(p+1) - t[,1]^(p+1))
  return(result)
}


fA_integral <- function(age, tb, p)
{
  retval <-1/(p+2)*(age^(p+2) - tb^(p+2)) -  tb/(p+1)*(age^(p+1) - tb^(p+1))
  return(retval)
}

lnSurv <- function(age, tb, te, b, p, beta1, alpha1, epsilon, psi, drate)
{
  acen <- matrix(gb.acen, nrow = gb.nRats, ncol = 1)
  tcen <- matrix(gb.tcen, nrow = gb.nRats, ncol = 1)
  dcen <- matrix(gb.dcen, nrow = gb.nRats, ncol = 1)
  lcen <- matrix(gb.lcen, nrow = gb.nRats, ncol = 1)
  
  time.matrix <- t
  time.matrix[,3] <- pmax(time.matrix[,3]- tlag, 1)

  k_constant <- beta1*exp(epsilon*(time.matrix[,3] - (t[,1] + (t[,2] - t[,1])/2) - tcen) + psi*(df$drate - dcen))
  
  lnS <- - 1/p*exp(b +p*log(time.matrix[,3]/lcen)) - exp(b)/p*k_constant*df$drate*(t[,2] - lcen)*exp(p*(t[,3]-lcen)) + exp(b)/(p^2)*k_constant*df$drate*exp(p*(time.matrix[,3]-lcen)) - exp(b)/(p^2)*k_constant*df$drate*exp(p*(t[,1]-lcen))
  
  return(exp(lnS))
}



# wrapper function for survival according to the ERR model
ERR_surv_vector <- function(t,par, drate, expdur)
{
  b <- par[1]
  p <- par[2]
  beta1 <- par[3]
  alpha1 <- par[4]
  epsilon <- par[5]
  psi <- par[6]
  tlag <- par[7]

    
  ndim <- length(t)

  
  m.t <- matrix(0, nrow= length(t), ncol = 3) # create internal time step matrix
  t <- pmax(t- tlag, 1)
  m.t[, 1] <- pmin(13, t-2/7)
  m.t[, 2] <- pmin(13 + expdur, t-1/7)
  m.t[, 3] <- t
  tsme <- m.t[,3] - (m.t[,1] + m.t[,2])/2
   
  lnS <- unlist(lapply(1:ndim, function(i) 
    lnSurv(m.t[i,3],m.t[i,1],m.t[i,2],b, p, beta1, alpha1, epsilon, psi, drate)))

    return (exp(lnS))
}


#----------------
# Define global exposure classes of Table 2
 
df$clsEER <- "unexposed"
df$clsEER[df$clsWLM4 == "low" & df$clsWL3 == "low"] <- "Ll"
df$clsEER[df$clsWLM4 == "medium" & df$clsWL3 == "low"] <- "Ml"
df$clsEER[df$clsWLM4 == "medium" & df$clsWL3 == "high"] <- "Mh"
df$clsEER[df$clsWLM4 == "high" & df$clsWL3 == "low"] <- "Hl"
df$clsEER[df$clsWLM4 == "high" & df$clsWL3 == "high"] <- "Hh"
df$clsEER <- fct_relevel(df$clsEER, "unexposed","Ll","Ml","Mh","Hl","Hh")
table(df$clsEER)

# Kaplan Meier plots
sv.ad <- survfit(Surv(aeofup/scfac,adeno) ~ clsEER, data = df)
sv.sq <- survfit(Surv(aeofup/scfac,epid) ~ clsEER, data = df)
sv.al <- survfit(Surv(aeofup/scfac,tumor) ~ clsEER, data = df)

# dataframe for ad
strata <- sv.ad$strata
names(strata) <- c("unexposed", "Ll", "Ml", "Mh", "Hl", "Hh")

newsurv <- sv.ad$surv
newnames <- c(rep(names(strata),strata))
ages <- sv.ad$time
km.ad <- data.frame(strata=newnames, age = ages, surv=newsurv)
km.ad$histo <- "ADC"
head(km.ad)

# dataframe for sq
strata <- sv.sq$strata
names(strata) <- c("unexposed", "Ll", "Ml", "Mh", "Hl", "Hh")

newsurv <- sv.sq$surv
newnames <- c(rep(names(strata),strata))
ages <- sv.sq$time
km.sq <- data.frame(strata=newnames, age = ages, surv=newsurv)
km.sq$histo <- "SQCC"
head(km.sq)



# dataframe for all

strata <- sv.al$strata
names(strata) <- c("unexposed", "Ll", "Ml", "Mh", "Hl", "Hh")

newsurv <- sv.al$surv
newnames <- c(rep(names(strata),strata))
ages <- sv.al$time
km.al <- data.frame(strata=newnames, age = ages, surv=newsurv)
km.al$histo <- "ALL LC"
head(km.al)

# bind all three together

pf.km <- rbind(km.al,km.ad,km.sq)
names(pf.km)[1] <- "expgrp"
names(pf.km)


pf.km$histo <- fct_relevel(pf.km$histo, "All LC","ADC","SQCC")
pf.km$expgrp <- fct_relevel(pf.km$expgrp, "unexposed","Ll","Ml","Mh","Hl","Hh")

head(pf.km)
str(pf.km)
summary(pf.km)

# predicted survival for adeno
ages <- pf.km$age
S.0 <- ERR_surv_vector(ages, wparadenoERR, 0, 0)
S.Ll <- ERR_surv_vector(ages, wparadenoERR, 70, 0.6)
S.Ml <- ERR_surv_vector(ages, wparadenoERR, 52, 8)
S.Mh <- ERR_surv_vector(ages, wparadenoERR, 565, 1)
S.Hl <- ERR_surv_vector(ages, wparadenoERR, 60, 51)
S.Hh <- ERR_surv_vector(ages, wparadenoERR, 390, 14)

hDd.ad <- data.frame(
  x = rep(ages, 1),
  y = c(S.0, S.Ll, S.Ml, S.Mh, S.Hl, S.Hh),
  expgrp = rep( c("unexposed","Ll","Ml","Mh","Hl","Hh" ), each = length(ages))
)

hDd.ad$histo <- "ADC"

# predicted survival for sq
ages <- pf.km$age
S.0 <- ERR_surv_vector(ages, wparsqERR, 0, 0)
S.Ll <- ERR_surv_vector(ages, wparsqERR, 70, 0.6)
S.Ml <- ERR_surv_vector(ages, wparsqERR, 52, 8)
S.Mh <- ERR_surv_vector(ages, wparsqERR, 565, 1)
S.Hl <- ERR_surv_vector(ages, wparsqERR, 60, 51)
S.Hh <- ERR_surv_vector(ages, wparsqERR, 390, 14)



hDd.sq <- data.frame(
  x = rep(ages, 1),
  y = c(S.0, S.Ll, S.Ml, S.Mh, S.Hl, S.Hh),
  expgrp = rep( c("unexposed","Ll","Ml","Mh","Hl","Hh" ), each = length(ages))
)

hDd.sq$histo <- "SQCC"


# predicted survival for all
ages <- pf.km$age
S.0 <- ERR_surv_vector(ages, wparallERR, 0, 0)
S.Ll <- ERR_surv_vector(ages, wparallERR, 70, 0.6)
S.Ml <- ERR_surv_vector(ages, wparallERR, 52, 8)
S.Mh <- ERR_surv_vector(ages, wparallERR, 565, 1)
S.Hl <- ERR_surv_vector(ages, wparallERR, 60, 51)
S.Hh <- ERR_surv_vector(ages, wparallERR, 390, 14)



hDd.al <- data.frame(
  x = rep(ages, 1),
  y = c(S.0, S.Ll, S.Ml, S.Mh, S.Hl, S.Hh),
  expgrp = rep( c("unexposed","Ll","Ml","Mh","Hl","Hh" ), each = length(ages))
)

hDd.al$histo <- "ALL LC"

# bind all three together

pf.md <- rbind(hDd.al,hDd.ad,hDd.sq)
names(pf.md)[1] <- "age"
names(pf.md)[2] <- "surv"
names(pf.md)


pf.md$histo <- fct_relevel(pf.md$histo, "ALL LC","ADC","SQCC")
pf.md$expgrp <- fct_relevel(pf.md$expgrp, "unexposed","Ll","Ml","Mh","Hl","Hh")

head(pf.md)
str(pf.md)
summary(pf.md)


fp.1 <- ggplot() + 
  geom_point(data = pf.km, aes(x=age, y=surv, color = expgrp, group = expgrp), size = 1) +
  geom_line(data = pf.md, aes(x=age, y=surv, color = expgrp, group = expgrp), linewidth = 0.5) +
  facet_grid(. ~ factor(histo, levels = c("ALL LC", "ADC", "SQCC"))) +
  scale_x_continuous(name = "Rat age (weeks)") +
  scale_y_continuous(name="Surviving fraction", limits=c(0,1), breaks = seq(0,1,0.2)) +
  scale_color_manual(values=cbPalette, name = "Exp. class") +
  theme(text = element_text(size=15),legend.position = c(.07,0.2)) 
print(fp.1)




#------------------------------------------------------------------------
# Excess relative risk function according to ERR model needed for Fig S7
#------------------------------------------------------------------------
hazard <- function(age, tb, te, b, p, beta1, alpha1, epsilon, psi, drate)
{
  acen <- matrix(gb.acen, nrow = gb.nRats, ncol = 1)
  tcen <- matrix(gb.tcen, nrow = gb.nRats, ncol = 1)
  dcen <- matrix(gb.dcen, nrow = gb.nRats, ncol = 1)
  lcen <- matrix(gb.lcen, nrow = gb.nRats, ncol = 1)
  
  time.matrix <- t
  time.matrix[,3] <- pmax(time.matrix[,3]- tlag, 1)

  k_constant <- beta1*exp(epsilon*(time.matrix[,3] - (t[,1] + (t[,2] - t[,1])/2) - tcen) + psi*(df$drate - dcen))
  
  haz <- exp(b + p*log(t[,3]/lcen))*(1 + k_constant*df$drate*(t[,2] - t[,1]))
  
  return(haz)
}



# wrapper function for ERR function according to the ERR model
ERR_vector <- function(t,par, drate, expdur)
{
  b <- par[1]
  p <- par[2]
  beta1 <- par[3]
  alpha1 <- par[4]
  epsilon <- par[5]
  psi <- par[6]
  tlag <- par[7]    
  ndim <- length(t)

  m.t <- matrix(0, nrow= length(t), ncol = 3) # create internal time step matrix
  t <- pmax(t- tlag, 1)
  m.t[, 1] <- pmin(13, t-2/7)
  m.t[, 2] <- pmin(13 + expdur, t-1/7)
  m.t[, 3] <- t
  tsme <- m.t[,3] - (m.t[,1] + m.t[,2])/2
  
  dose <- drate*expdur 
  haz0 <- unlist(lapply(1:ndim, function(i) 
    hazard(    lnSurv(m.t[i,3],m.t[i,1],m.t[i,2],b, p, beta1, alpha1, epsilon, psi, 0)))
  haz1 <- unlist(lapply(1:ndim, function(i) 
    hazard(    lnSurv(m.t[i,3],m.t[i,1],m.t[i,2],b, p, beta1, alpha1, epsilon, psi, drate)))
  
  retval <- haz/haz0 -1
  return(retval/dose)

}

tlag <- wpar[7]

p <- ggplot(data.frame(x = c(5,175)), aes(x = x))+ 
  geom_path(aes(colour= cbPalette[2]), stat="function", fun = ERR_vector, args = list(par =wpar, drate = 390, expdur = 14), size = 1.1)+
  #geom_path(aes(colour= cbPalette[3]), stat="function", fun = ERR_ERR_vector, args = list(par =wparadenoERR, drate = 10, expdur = 50), size = 1.1)+
  geom_vline(xintercept = 13, linetype = 2,  size = 0.9) +
  annotate(geom = "text", x = 9, y = 10/500, label = "Begin of Exposure",angle = 90)+
  
  geom_vline(xintercept = 27, linetype = 2, col = cbPalette[2], size = 0.9) +
  annotate(geom = "text", x = 27, y = 10/500, label = "End of Exposure", col = cbPalette[2], angle = 90)+
  
  geom_vline(xintercept = 27 + tlag, linetype = 2, col = cbPalette[2], size = 0.9) +
  geom_segment(aes(x = 27, y = 13.5/500, xend = 27+tlag, yend = 13.5/500), col = cbPalette[2], size= 0.9,
               arrow = arrow(length = unit(0.2, "cm"), ends = "both")) +
    xlab("Rat Age (weeks)") +
  ylab("Excess Relatie Risk/WLM") +
coord_cartesian(xlim = c(5, 175), ylim = c(0, 1), expand = FALSE)+
  theme(legend.text = element_text(
    #face = "bold",
    size = 14)
  ) +
  theme(axis.title.x = element_text( size = 14)) +
  #ggtitle ("ADC")+
  theme(axis.title.y = element_text( size = 14),legend.position=c(.8,.9))

p


