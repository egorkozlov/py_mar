#get the final sample for the analysis


## remove (almost) everything in the working environment.
rm(list = ls())


my_packages <- c("survival","rms","plm","lmtest","stargazer","foreing","nnet","ggplot2","reshape2","haven", "mlogit", "margins", "gmnl","VGAM","MNP","sampleSelection","")                                        # Specify your packages
not_installed <- my_packages[!(my_packages %in% installed.packages()[ , "Package"])]    # Extract not installed packages
if(length(not_installed)) install.packages(not_installed)                               # Install not installed packages


#for reading .dta file
library(haven)

#for multinomial logit+tobit+probit
library(mlogit)
library(mnlogit)
library(margins)
library(gmnl)
library(VGAM)
library(MNP)
library(sampleSelection)
library(survival)
library(rms)

library(multiwayvcov)
library(margins)

#For testing if parameters are different
library(car)

#fixed effects regression
library(plm)
require(lmtest)

#for exporting in Latex
library(stargazer)
require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)

library(did)
library(DRDID)
##################################################################
##################################################################
####  SAMPLE SELECTION AS A FIRST THING
##################################################################
##################################################################
directory<-"C:/Users/Fabio/Dropbox/JMP/empirical analysis/Analysis/duration"
#directory<-"C:\\Users\\Fabio\\Dropbox\\cohabitation and divorce laws\\coh_text"
setwd(directory)

#Here Sample selection: note that some specific observation
#might be dropped for different reasons depending on the
#Analysis we do
select<-function(data,age_min,age_max,femalen) {
  
  data<-data[data$age>=age_min,]
  data<-data[data$age<=age_max,]
  data<-data[data$female>=femalen,]
  
}

##################################################################
##################################################################
####  MARRIAGE versus COHABITATION
##################################################################
##################################################################
#Rules of thumb for comparing multinomial logit and multinomial probit coefficients
ind <-read_dta("C:/Users/Fabio/Dropbox/JMP/empirical analysis/Analysis/duration/NSFG88.dta")

G <- length(unique(ind$stat))
c <- G/(G - 1)

ind$failm<-0
ind$failm[ind$fail2==1]<-1

ind$failc<-0
ind$failc[ind$fail2==2]<-1

lm1 <- lm(failm ~ unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          data=ind,#weight=wgt_n,
          na.action=na.exclude)


lm2 <- lm(failm ~  unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          data=ind[ind$keep==1,],#weight=wgt_n,
          na.action=na.exclude)

lm3 <- lm(failm ~  unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          index=c('state'),
          data=ind[ind$nsfh==1,],#weight=wgt_n,
          na.action=na.exclude)

lm4 <- lm(failm ~  unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          index=c('state'), 
          data=ind[ind$nsfh==0,],#weight=wgt_n,
          na.action=na.exclude)

#put the stuff in latex
v1<-multiwayvcov::cluster.vcov(lm1, 
                               cbind(ind$stat, ind$date_rel),
                               use_white = F, 
                               df_correction = F)

v2<-multiwayvcov::cluster.vcov(lm2, 
                               cbind(ind[ind$keep==1,]$stat, ind[ind$keep==1,]$date_rel),
                               use_white = F, 
                               df_correction = F)

v3<-multiwayvcov::cluster.vcov(lm3, 
                               cbind(ind[ind$nsfh==1,]$stat, ind[ind$nsfh==1,]$date_rel),
                               use_white = F, 
                               df_correction = F)

v4<-multiwayvcov::cluster.vcov(lm4, 
                               cbind(ind[ind$nsfh==0,]$stat, ind[ind$nsfh==0,]$date_rel),
                               use_white = F, 
                               df_correction = F)

Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=FALSE,
              se = list(coeftest(lm1,v1)[, 2],
                        coeftest(lm2,v2)[, 2],
                        coeftest(lm3,v3)[, 2],
                        coeftest(lm4,v4)[, 2]),
              omit.table.layout = "n",
              dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
              covariate.labels=c("Unilateral Divorce"),
              dep.var.labels=c("","","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              keep.stat=c("n","rsq"),no.space=TRUE,
              keep = c('unid'),
              font.size='footnotesize',
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Duration Polynomial", "Yes", "Yes", "Yes", "Yes"),
                               c("Socio-economic controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\linm.tex',sep=""))



lm1 <- lm(failc ~ unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          data=ind,#weight=wgt_n,
          na.action=na.exclude)


lm2 <- lm(failc ~  unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          data=ind[ind$keep==1,],#weight=wgt_n,
          na.action=na.exclude)

lm3 <- lm(failc ~  unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          index=c('state'),
          data=ind[ind$nsfh==1,],#weight=wgt_n,
          na.action=na.exclude)

lm4 <- lm(failc ~  unid+factor(state)+factor(beg)+t+t2+t3+coll+female,
          index=c('state'), 
          data=ind[ind$nsfh==0,],#weight=wgt_n,
          na.action=na.exclude)

#put the stuff in latex
v1<-multiwayvcov::cluster.vcov(lm1, 
                               cbind(ind$stat, ind$date_rel),
                               use_white = F, 
                               df_correction = F)

v2<-multiwayvcov::cluster.vcov(lm2, 
                               cbind(ind[ind$keep==1,]$stat, ind[ind$keep==1,]$date_rel),
                               use_white = F, 
                               df_correction = F)

v3<-multiwayvcov::cluster.vcov(lm3, 
                               cbind(ind[ind$nsfh==1,]$stat, ind[ind$nsfh==1,]$date_rel),
                               use_white = F, 
                               df_correction = F)

v4<-multiwayvcov::cluster.vcov(lm4, 
                               cbind(ind[ind$nsfh==0,]$stat, ind[ind$nsfh==0,]$date_rel),
                               use_white = F, 
                               df_correction = F)

Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=FALSE,
              se = list(coeftest(lm1,v1)[, 2],
                        coeftest(lm2,v2)[, 2],
                        coeftest(lm3,v3)[, 2],
                        coeftest(lm4,v4)[, 2]),
              omit.table.layout = "n",
              dep.var.caption=c("\\textit{Dependent variable: Broke up (0/1)}"),
              covariate.labels=c("Unilateral Divorce"),
              dep.var.labels=c("","","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              keep.stat=c("n","rsq"),no.space=TRUE,
              keep = c('unid'),
              font.size='footnotesize',
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Duration Polynomial", "Yes", "Yes", "Yes", "Yes"),
                               c("Socio-economic controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\linc.tex',sep=""))
