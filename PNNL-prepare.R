#----------------------------------------------------------
# Prepare data set of PNNL rats
# jck, 2022/09/08
#----------------------------------------------------------
rm(list = ls()) # remove all objects from the current workspace

#library(msce) # Cristoforo's package
library(forcats)
library(ggplot2)

#------------------------------------------------------------------
# read raw data newly filtered from rat.wh
# from Heidenreich et al. RasRes 151, 209-17, 1999
#------------------------------------------------------------------
setwd("I:/Rats/stats") # project root directory
setwd("I:/Rats/stats/data")
lungCancerRadon <- read.csv(file = "ratwh.csv")
setwd("I:/Rats/stats") 

is.data.frame(lungCancerRadon)
df <- lungCancerRadon
head(df)
dim(df)[1] # 4260

# rescale WLM, WL
df$hWLM <- df$dose/100
df$hWL <- df$drate/100

# binary
df$exposed <- "no"
df$exposed[df$dose > 0] <- "yes"

# adjust deathcode
table(df$deathcode)

df$deathcode[df$deathcode == " A"] <- "AK"
df$deathcode[df$deathcode == " K"] <- "AK"
df$deathcode[df$deathcode == " S"] <- "S"
df$deathcode[df$deathcode == " E"] <- "E"
df$deathcode[df$deathcode == " D"] <- "D"

table(df$deathcode)
df$deathcode <- fct_relevel(df$deathcode,"D","S","E","AK")
table(df$deathcode)

#-------------------------------------------------
# Exposure classes
#-------------------------------------------------
# define cumulative exposure class (WLM)
df$clsWLM <- "unexposed"
df$clsWLM[df$dose > 0] <- "20"
df$clsWLM[df$dose > 30] <- "40"
df$clsWLM[df$dose > 60] <- "80"
df$clsWLM[df$dose > 120] <- "160"
df$clsWLM[df$dose > 300] <- "320"
df$clsWLM[df$dose > 600] <- "640"
df$clsWLM[df$dose > 1000] <- "1280"
df$clsWLM[df$dose > 2000] <- "2560"
df$clsWLM[df$dose > 5000] <- "5120"
df$clsWLM[df$dose > 10000] <- "10240"
df$clsWLM <- fct_relevel(df$clsWLM, "unexposed", "20", "40","80","160","320","640","1280","2560","5120","10240")
table(df$clsWLM)

# define cumulative exposure class (WLM), 4 levels
df$clsWLM4 <- "unexposed"
df$clsWLM4[df$dose > 0] <- "low"
df$clsWLM4[df$dose > 100] <- "medium"
df$clsWLM4[df$dose > 1000] <- "high"
df$clsWLM4 <- fct_relevel(df$clsWLM4, "unexposed", "low", "medium", "high")
table(df$clsWLM4)

# derived from actual dose rate
table(df$drate)
df$clsWL <- "unexposed"
df$clsWL[df$drate > 0] <- "5"
df$clsWL[df$drate  > 25] <- "50"
df$clsWL[df$drate  > 100] <- "200"
df$clsWL[df$drate  > 400] <- "500"
df$clsWL[df$drate  > 650] <- "750"
df$clsWL <- fct_relevel(df$clsWL, "unexposed", "5", "50","200","500","750")
table(df$clsWL)
aggregate(df$drate,list(df$clsWL),mean)
aggregate(df$drate,list(df$clsWL),min)
aggregate(df$drate,list(df$clsWL),max)

# derived from actual dose rate, three levels
table(df$drate)
df$clsWL3 <- "unexposed"
df$clsWL3[df$drate > 0] <- "low"
df$clsWL3[df$drate  > 100] <- "high"
df$clsWL3 <- fct_relevel(df$clsWL3, "unexposed", "light", "heavy")
table(df$clsWL3)

table(df$clsWL,df$clsWL3)
#-------------------------------------------------
# Fatal tumors
#-------------------------------------------------
table(df$tumor)
df$tumfat <- 0
df$tumfat[df$TMR_ADENCAR > 2] <- 1
df$tumfat[df$TMR_EPID > 2] <- 1
df$tumfat[df$TMR_ADENOSQ > 2] <- 1
df$tumfat[df$TMR_SARCOMA > 2] <- 1
table(df$tumfat)

#-------------------------------------------------
# Tumor histology
#-------------------------------------------------
df$adeno <- 0
df$adeno[df$TMR_ADENCAR > 0] <- 1
df$epid <- 0
df$epid[df$TMR_EPID > 0] <- 1
df$adesq <- 0
df$adesq[df$TMR_ADENOSQ > 0] <- 1
df$sarco <- 0
df$sarco[df$TMR_SARCOMA > 0] <- 1

sum(df$adeno) # 302
sum(df$epid) # 104
sum(df$adesq) # 36
sum(df$sarco) # 35

#-------------------------------------------------
# Age groups (months)
#-------------------------------------------------
dpm <- 365/12 # days per month
df$agrp <- "<5"
df$agrp[df$aeofup/dpm > 4] <- "5-8"
df$agrp[df$aeofup/dpm > 8] <- "9-12"
df$agrp[df$aeofup/dpm > 12] <- "13-16"
df$agrp[df$aeofup/dpm > 16] <- "17-20"
df$agrp[df$aeofup/dpm > 20] <- "21-24"
df$agrp[df$aeofup/dpm > 24] <- "25-28"
df$agrp[df$aeofup/dpm > 28] <- "29-32"
df$agrp[df$aeofup/dpm > 32] <- "33-36"
df$agrp[df$aeofup/dpm > 36] <- "37-40"
df$agrp[df$aeofup/dpm > 40] <- ">40"
df$agrp <- fct_relevel(df$agrp, "<5", "5-8", "9-12","13-16","17-20","21-24","25-28","29-32","33-36","37-40",">40")
table(df$agrp)

aggregate(df$tumor,list(df$agrp),sum)
aggregate(df$tumor,list(df$clsWLM),sum)
aggregate(df$tumor,list(df$clsWLM),length)
round(aggregate(df$tumor,list(df$clsWLM),sum)$x/aggregate(df$tumor,list(df$clsWLM),length)$x*100,1)


#----------------------------------------
# correct for death during exposure
#----------------------------------------
df$dae <- FALSE # death after exposure
df$dde <- TRUE # death during exposure
df$dae[(df$aeofup - df$aeoexp) > 0] <- TRUE
df$dde[(df$aeofup - df$aeoexp) > 0] <- FALSE
table(df$dae,df$dde)
#       FALSE TRUE
# FALSE     0   53
# TRUE   4207    0

# adjust ages to meet msce matrix format
df$aeofup[df$dde == T] <- df$aeoexp[df$dde == T]
df$aeoexp[df$dde == T] <- df$aboexp[df$dde == T]
df$aboexp[df$dde == T] <- 0

df$aboexp[df$dde == T]
df$aeoexp[df$dde == T]
df$aeofup[df$dde == T]

#-------------------------------------------------
# Assemble data frame for saving
#-------------------------------------------------
names(df)
df.pnnl <- df[,c("id", "aboexp","aeoexp","aeofup","expdur","agrp","exposed",
                 "dose","drate","hWLM", "hWL","clsWLM","clsWLM4","clsWL","clsWL3",
                 "tumor","tumfat","adeno","epid","adesq","sarco","deathcode","dae","dde")]
names(df.pnnl)
dim(df.pnnl)
str(df.pnnl)
summary(df.pnnl)

setwd("I:/Rats/stats/data")
save(df.pnnl,file = "PNNLrats-231121.Rdata")

setwd("I:/Rats/stats")
