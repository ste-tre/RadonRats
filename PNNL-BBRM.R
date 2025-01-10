# Fit BBRM model for PNNL rats
# st 01/02/24
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

#---------------------------------------------------------------------------------
# set data frame for analysis
#---------------------------------------------------------------------------------

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

# msce_numerical
gb.nInnerSteps <- 1000 # default

nInnerSteps <- 1000 # inner steps for numerical fitting
#-----------------------------------------------------------------------------
# select models and likelihoods
#-----------------------------------------------------------------------------
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

#-----------------------------------------------------------------------------
# import likelihood functions
#-----------------------------------------------------------------------------
source("I:/Rats/stats/subscripts/likelihood.R")

# controls CI printing
has.profile <- FALSE
has.Hessian <- FALSE

smloglik <- msce3_mle2_loglik_IP
snloglik <- msce3_nlm_loglik_IP
sdeviance <- msce3_deviance_IP

#---------------------------
# read estimated parameters
#---------------------------

parDF <- read.csv("I:/Rats/stats/results/tsce+TG/Best_parms/tsce-TG-all-all-bio_parms.csv")

setwd("I:/Rats/stats")
upar <- vector()
upar <- as.numeric(parDF$parval[c(1:8)])
upar

wpar <- vector() # working parameters, will be modified
wpar[1:8] <- upar[1:8]

names(wpar) <- c("X0", "Y1", "Y2", "gam0", "gam1", "gam2", "alpha1", "gamc")
wpar
#-------------------------------------------------------------------------
# fit parameters with bbmle package 

sdeviance(wpar[1], wpar[2], wpar[3], wpar[4], wpar[5], wpar[6], wpar[7], wpar[8], df) 

mmcall = 0
durat<- system.time(
  mle.3IP <- mle2(minuslog = smloglik,
                  start=list(X0 = wpar[1], Y1 = wpar[2], Y2= wpar[3], gam0 = wpar[4], 
                             gam1 = wpar[5], gam2 = wpar[6], alpha1 = wpar[7], gamc = wpar[8]), 
                  parameters=list(X0~1, Y1~1, Y2~2, gam0~1, gam1~1, gam2~1, alpha1~1, gamc ~1), 
                  fixed=list(
                    #X0 = wpar[1],
                    #Y1 = wpar[2],
                    Y2 = 0 #wpar[3]
                    #gam0 = wpar[4]
                    #gam1 = wpar[5],
                    #gam2 = wpar[6]
                    #alpha1 = wpar[7]
                    #gamc = wpar[8]
                  ),
                  
                  lower = c(X0 = 1.5, Y1 = 0.1,  gam0 = 0.03, gam1 = 5e-2, gam2 = 0.01, alpha1 = 1,  gamc = 0.01),
                  upper = c(X0 = 2.5, Y1 = 1, gam0 = 0.05, gam1 = 5e-1, gam2 = 0.1, alpha1 = 10, gamc = 0.3),
                  
                  method = "L-BFGS-B",
                  #method = "CG",
                  #method = "Brent",
                  data = df)
) # system.time
durat

cat(sprintf("Total calls: %d", mmcall))

summary(mle.3IP)

wpar <- coef(mle.3IP)
names(wpar) <- c("X0", "Y1", "Y2", "gam0", "gam1", "gam2", "alpha1", "gamc")

#-----------------------------------------------------------------------------------
# Hazard, Survival and ERR/WLM functions used in Figures 3, 6, S2, S7
#-----------------------------------------------------------------------------------
# Hazard as function of age and drate
msce_hazard <- function(t, par, drate, tb = 13,  expdur){
  
  # time information in weeks!
  alpha0 <- matrix(1, nrow = length(t), ncol = 3) 
  
  X0 <- exp(par[1])
  Y1 <- exp(par[2])
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  alpha1 <- par[7]
  gamc <- par[8]
  
  del0 <- exp(-12)
  M <- exp(-13.82) 
  alphac <- matrix(alpha1, nrow = length(t), ncol = 3)
  
  before <- t <= 13
  during <- 13 < t & t <= 13 + expdur
  after <- t > 13 + expdur
  time_matrix <- matrix(0, nrow = length(t), ncol = 3)
  
  time_matrix[, 1] <- 13*after
  time_matrix[, 2] <- 13*during + after*(13 + expdur)
  time_matrix[, 3] <- pmax(t, 1)

  
  expInd <- matrix(0, nrow = length(t), ncol = 3)
  expInd[, 2] <- after
  expInd[, 3] <- during  

  # biological rates
  Nnu0 <- X0*alpha0 + X0*alpha0*Y1*drate*exp(-Y2*drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=length(t),ncol=3)
  # detection rate
  nu2 <- matrix(M, nrow = length(t), ncol = 3) # M/alpha1
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow= length(t),ncol=3) + gam2*(1-exp(-gam1/gam2*drate))*expInd
  gammac <- matrix(gamc, nrow = length(t), ncol = 3)
  
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha0,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  
  result <- msce_numerical(time_matrix, parList, innerSteps = 1000)
  haz <- result$hazard
  return(haz)
}
#-------------------------------------------------------------------
# Survival function
msce_survival <- function(t, par, drate, tb = 13, expdur){
  
  # time information in weeks!
  alpha0 <- matrix(1, nrow = length(t), ncol = 3)
  
  X0 <- exp(par[1])
  Y1 <- exp(par[2])
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  alpha1 <- par[7]
  gamc <- par[8]
  
  del0 <- exp(-12)
  M <- exp(-13.82) 
  alphac <- matrix(alpha1, nrow = length(t), ncol = 3)
  
  before <- t <= 13
  during <- 13 < t & t <= 13 + expdur
  after <- t > 13 + expdur
  time_matrix <- matrix(0, nrow = length(t), ncol = 3)
  
  time_matrix[, 1] <- 13*after
  time_matrix[, 2] <- 13*during + after*(13 + expdur)
  time_matrix[, 3] <- pmax(t, 1)

  
  expInd <- matrix(0, nrow = length(t), ncol = 3)
  expInd[, 2] <- after
  expInd[, 3] <- during  
  
  # biological rates
  Nnu0 <- X0*alpha0 + X0*alpha0*Y1*drate*exp(-Y2*drate)*expInd
  nu1 <- matrix(del0/alpha0,nrow=length(t),ncol=3)
  # detection rate
  nu2 <- matrix(M, nrow = length(t), ncol = 3) # M/alpha1
  # Net growth intermediate and malignant cells
  gamma <- matrix(gam0,nrow= length(t),ncol=3) + gam2*(1-exp(-gam1/gam2*drate))*expInd
  gammac <- matrix(gamc, nrow = length(t), ncol = 3)
  
  parList = list(Nnu0=Nnu0,nu1=nu1,alpha1=alpha0,gamma1=gamma, alpha2 = alphac, gamma2 = gammac, nu2=nu2) 
  
  result <- msce_numerical(time_matrix, parList, innerSteps = 1000)
  lnS <- result$lnSurvival
  return(exp(lnS))
}

#------------------------------------------
# ERR per WLM function
msce_ERR <- function(t, par, drate, tb = 13,  expdur){
	h0 <- msce_hazard(t, par, drate = 0, tb = 0,  expdur = 0) # hazard without exposure
	h1 <- msce_hazard(t, par, drate, tb = 13,  expdur)
	exposure_total <- drate*expdur/7.08
	retval <- h1/h0 -1
	return(retval/exposure_total)
}

#----------------------------------------------------------------
# Comparison of Kaplan-Meier plots and survival curves for Fig 3
#----------------------------------------------------------------
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
S.0 <- msce_survival(ages, wparadeno, 0, 0)
S.Ll <- msce_survival(ages, wparadeno, 70, 0.6)
S.Ml <- msce_survival(ages, wparadeno, 52, 8)
S.Mh <- msce_survival(ages, wparadeno, 565, 1)
S.Hl <- msce_survival(ages, wparadeno, 60, 51)
S.Hh <- msce_survival(ages, wparadeno, 390, 14)

hDd.ad <- data.frame(
  x = rep(ages, 1),
  y = c(S.0, S.Ll, S.Ml, S.Mh, S.Hl, S.Hh),
  expgrp = rep( c("unexposed","Ll","Ml","Mh","Hl","Hh" ), each = length(ages))
)

hDd.ad$histo <- "ADC"

# predicted survival for sq
ages <- pf.km$age
S.0 <- msce_survival(ages, wparsq, 0, 0)
S.Ll <- msce_survival(ages, wparsq, 70, 0.6)
S.Ml <- msce_survival(ages, wparsq, 52, 8)
S.Mh <- msce_survival(ages, wparsq, 565, 1)
S.Hl <- msce_survival(ages, wparsq, 60, 51)
S.Hh <- msce_survival(ages, wparsq, 390, 14)

hDd.sq <- data.frame(
  x = rep(ages, 1),
  y = c(S.0, S.Ll, S.Ml, S.Mh, S.Hl, S.Hh),
  expgrp = rep( c("unexposed","Ll","Ml","Mh","Hl","Hh" ), each = length(ages))
)

hDd.sq$histo <- "SQCC"


# predicted survival for all
ages <- pf.km$age
S.0 <- msce_survival(ages, wparall, 0, 0)
S.Ll <- msce_survival(ages, wparall, 70, 0.6)
S.Ml <- msce_survival(ages, wparall, 52, 8)
S.Mh <- msce_survival(ages, wparall, 565, 1)
S.Hl <- msce_survival(ages, wparall, 60, 51)
S.Hh <- msce_survival(ages, wparall, 390, 14)

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
  theme(text = element_text(size=15),legend.position = c(.1,0.2)) 
print(fp.1)

#-------------------------------------------
# Figure 6
#-------------------------------------------

ages <- seq(0, 175, by = 1)

E0 <- msce_hazard(ages, wparsq, drate = 0, expdur = 0 )
E.Ml <- (msce_hazard(ages, wparsq, drate = 52, expdur = 416/52)/E0 -1)/(416)
E.Mh <- (msce_hazard(ages, wparsq, drate = 565, expdur = 1)/E0-1)/(565)

Sq.Ml <- data.frame(
  ages = ages,
  y = E.Ml,
  exp = rep("Medium-light", length(ages)),
endexp = 21, #64,
max = 63,
y1 = 0.001,
x1 = 17,
x2 = (21 + 63)/2,
y2 = 0.002, 
y3 =  1.05,
y4 = 0.002,
TC = 53,
TD = 110,
  histo = rep("SQCC", length(ages))

)

Sq.Mh <- data.frame(
  ages = rep(ages, 1),
  y = E.Mh,
  exp = rep("Medium-heavy", length(ages)),
endexp = 14, 
max = 56,
y1 = 0.001,
x1 = 10,
x2 = (25+ 56)/2,
y2 = 0.002,
y3 = 1.05,
y4 = 0.002,
TC = 48,
TD = 114,
  histo = rep("SQCC", length(ages))
)

# adeno
ages <- seq(0,175, by = 1)

E0 <- msce_hazard(ages, wparadeno, drate = 0, expdur = 0 )
E.Ml <- (msce_hazard(ages, wparadeno, drate = 52, expdur = 8)/E0 -1)/(416)
E.Mh <- (msce_hazard(ages, wparadeno, drate = 565, expdur = 1)/E0-1)/(565)

Ad.Ml <- data.frame(
  ages = ages,
  y = E.Ml,
  exp = rep("Medium-light", length(ages)),
endexp = 21,
max = 75,
y1 = 0.04,
x1 = 17, 
x2 = (21+75)/2,
y2 = 0.06,
y3 = 1.05, 
y4 = 0.005,
TC = 34,
TD = 112,
  histo = rep("ADC", length(ages))

)

Ad.Mh <- data.frame(
  ages = rep(ages,1),
  y = E.Mh,
  exp = rep("Medium-heavy", length(ages)),
endexp = 14, 
max = 68,
y1 = 0.04,
x1 = 10,
x2 = (14+68)/2,
y2 = 0.06,
y3 = 1.05,
y4 = 0.005,
TC = 22,
TD = 104,
  histo = rep("ADC", length(ages))

)


pf.md <- rbind(Sq.Ml, Sq.Mh, Ad.Ml, Ad.Mh)
names(pf.md)

pf.md$histo <- fct_relevel(pf.md$histo, "ADC","SQCC")
pf.md$exp <- fct_relevel(pf.md$exp, "Medium-light","Medium-heavy" )

pf.md$cluster <- paste(pf.md$exp, pf.md$histo, sep = "-")

pf.md$cluster <- fct_relevel(pf.md$cluster,  "Medium-light-SQCC", "Medium-heavy-SQCC", "Medium-light-ADC", "Medium-heavy-ADC" )

head(pf.md)
str(pf.md)
summary(pf.md)

library(latex2exp)

fp.1 <- ggplot() + 
  geom_line(data = pf.md, aes(x=ages, y=y, col = exp), linewidth = 0.9) +
  facet_grid(histo ~ exp, scales = "free_y") +
  xlab("Rat age (weeks)") +
  ylab("ERR per WLM") +
  #coord_cartesian(xlim = c(0, 175), expand = FALSE) +
  geom_vline(data = pf.md, aes(xintercept = endexp, col = exp ) , linetype = 2)+
  geom_text(data = pf.md, aes(x = x1, y = y1, angle = 90), label = "End of exposure", check_overlap = T , size = 4)+
  geom_vline(data = pf.md, aes(xintercept = max, col = exp), linetype = 2)+
  geom_segment(data = pf.md, aes(x = endexp, y = y2, xend = max, yend = y2, col = exp),arrow = arrow(length = unit(0.2, "cm"), ends = "both") ) +
  geom_text( data = pf.md, aes(x = x2, y = y3*y2), label = TeX("Mean tumor age $A_T$"), check_overlap = T, size = 3.5)+
  scale_color_manual(values=cbPalette[3:4], name = "Exposure") +
  theme(text = element_text(size=15),legend.position = c(0.9,0.5)) 

print(fp.1)

#--------------------------------------------------------
# Figure S2: Net growth rate as function of exposre rate
#--------------------------------------------------------

gammas <- function(par, drate) {
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  
  gamma <- gam0 + gam2*(1-exp(-gam1/gam2*drate))
  return(gamma)
}

p <- ggplot(data.frame(x = c(min(df$drate), max(df$drate))), aes(x = x))+
  geom_path(aes(colour= cbPalette[3]), stat="function", fun = gammas, args = list(par = wparsqlog[1:6]), size = 1.1)+
  geom_path(aes(colour= cbPalette[4]), stat="function", fun = gammas, args = list(par = wparsq[1:6]), size = 1.1)+
  geom_path(aes(colour= cbPalette[2]), stat="function", fun = gammas, args = list(par = wparadeno[1:6]), size = 1.1)+
  geom_path(aes(colour= cbPalette[1]), stat="function", fun = gammas, args = list(par = wparall[1:6]), size = 1.1)+
  coord_cartesian(xlim = c(min(df$drate), max(df$drate)), ylim = c(0,0.5), expand = FALSE)+
  scale_colour_identity("Subgroup", guide="legend", 
                        labels = c( "SQCC (I+P)","SQCC", "ADC",  "All LC"), 
                        breaks = cbPalette[c(3, 4,2,1)]) +
  xlab("Exposure Rate (WL)") +
  ylab("Net Growth rate (cells/week)") +
  theme(legend.text = element_text(
    #face = "bold",
    size = 12),
    legend.position=c(.9,.45)
  ) +
  theme(axis.title.x = element_text( size = 12)) +
  theme(axis.title.y = element_text( size = 12))
p

#-------------------------------------------
# Expected and observed cases for figure S5
#-------------------------------------------

# msce applied to the dataset -> lnS for each rat
msce_rats <- function(par)
{
  mdl_call()
  
  # assume alpha to be small and constant
  alpha0 <- gb.alpha0
  alpha <- matrix(alpha0,nrow=gb.nRats,ncol=gb.nTimesteps)
  
  
  X0 <- exp(par[1])
  Y1 <- par[2]
  Y2 <- par[3]
  gam0 <- par[4]
  gam1 <- par[5]
  gam2 <- par[6]
  alpha1 <- par[7]
  gamc <- par[8]
  
  del0 <- exp(-12)#exp(-8.79)#10^-6
  M <-  exp(-13.82) #exp(-18.42) #10^-7 #exp(-14.4)#
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
  
  
  
  return (lnS)
}


lnS <- msce_rats(wpar)
sum(-lnS)
sum(df$tumor)

age_class <- numeric(length = 25)
expected <- numeric(length = 25)
observed <- numeric(length = 25)
ptui <- numeric(length = 25)
ptli <- numeric(length = 25)
age_name <- vector(mode="character", length=25)
for (iter in 1:24){
  age_class[iter] <- (iter-1)*50 +25
  age_name[iter] <- paste(as.character((iter-1)*50), as.character((iter)*50), sep = "-")
  expected[iter] <- sum(-lnS[df$aeofup >= (iter-1)*50 & df$aeofup < (iter)*50] ) # expected
  observed[iter] <- sum(df$tumor[df$aeofup >= (iter-1)*50 & df$aeofup < (iter)*50 ]) # observed
  
  ptui[iter] <- poisson.test(observed[iter], T = 1, alternative = "two.sided", conf.level = 0.95)$conf.int[2]
  ptli[iter] <- poisson.test(observed[iter], T = 1, alternative = "two.sided", conf.level = 0.95)$conf.int[1]
}
age_class[25] <- 1225
age_name[25] <- ">1200"
expected[25] <- sum(1-lnS[df$aeofup >= 1200] ) #expected
observed[25] <- sum(df$tumor[df$aeofup >= 1200 ]) # observed

ptui[25] <- poisson.test(observed[25], T = 1, alternative = "two.sided", conf.level = 0.95)$conf.int[2]
ptli[25] <- poisson.test(observed[25], T = 1, alternative = "two.sided", conf.level = 0.95)$conf.int[1]

-sum(lnS)

exp_obs <- data.frame(
  x = rep(age_class/7, 1),
  y = c(expected,observed),
  Age_name = age_name,
  exporobs = rep(c("Expected", "Observed"), each = 25)
)

p <- ggplot(data = exp_obs) +
  #geom_point(data = exp_obs[exp_obs$exporobs == "Expected",], aes(x = x, y= y), col = cbPalette[4], shape = 16, size = 2)+
  #geom_point(data = exp_obs[exp_obs$exporobs == "Observed",], aes(x = x, y= y), col = "black", shape = 4)+
  geom_point(aes(x = x, y = y, colour = factor(exporobs)))+
  scale_color_manual(name = "", breaks = c("Expected", "Observed"), values = cbPalette[c(2, 4)])+
  geom_line(data = exp_obs[exp_obs$exporobs == "Observed",], aes(x = x, y= y))+
  geom_errorbar(data = exp_obs[exp_obs$exporobs == "Observed",], aes(x = x, y= y, ymin= ptli, ymax= ptui), colour="black", width=5)+
  xlab("Age (weeks)") + ylab("Number of Lung carcinomas")

print(p)  


# exposure classes
exp_class <- character()
rate_class <- character()
y <- numeric()
w <- levels(df$clsWLM)
u <- levels(df$clsWL)
z <- numeric()
ptui <- numeric()
ptli <- numeric()
xlb <- character()
glob_class <- numeric()


index_y <- 0

for (iter in 1:length(w)){
  for (index in 1: length(u)) {
    if (length(df$id[df$clsWLM == w[iter] & df$clsWL == u[index]])){
      index_y <- index_y +1
      exp_class[index_y] <- w[iter]
      rate_class[index_y] <- u[index]
      glob_class[index_y] <- index_y
      y[index_y] <- sum(-lnS[df$clsWLM == w[iter] & df$clsWL == u[index]])
      z[index_y] <- sum(df$tumor[df$clsWLM == w[iter] & df$clsWL == u[index] ])
      ptui[index_y] <- poisson.test(z[index_y], T = 1, alternative = "two.sided", conf.level = 0.95)$conf.int[2]
      ptli[index_y] <- poisson.test(z[index_y], T = 1, alternative = "two.sided", conf.level = 0.95)$conf.int[1]
      xlb[index_y] <- paste(round(ptli[index_y]), round(ptui[index_y]), sep= "-")
      
    }
  }
  
  
  
  
}

expected <- y
observed <- z

exp_obs <- data.frame(
  x = rep(glob_class, 1),
  y = c(expected,observed),
  exporobs = rep(c("Expected", "Observed"), each = 22)
)

p <- ggplot(data = exp_obs) +
  #geom_point(data = exp_obs[exp_obs$exporobs == "Expected",], aes(x = x, y= y), col = cbPalette[4], shape = 16, size = 2)+
  #geom_point(data = exp_obs[exp_obs$exporobs == "Observed",], aes(x = x, y= y), col = "black", shape = 4)+
  geom_point(aes(x = x, y = y, colour = factor(exporobs)))+
  scale_color_manual(name = "", breaks = c("Expected", "Observed"), values = cbPalette[c(2, 4)])+
  geom_line(data = exp_obs[exp_obs$exporobs == "Observed",], aes(x = x, y= y))+
  geom_errorbar(data = exp_obs[exp_obs$exporobs == "Observed",], aes(x = x, y= y, ymin= ptli, ymax= ptui), colour="black", width=0.5)+
  xlab("Exposure class") + ylab("Number of Lung carcinomas")

print(p)  

