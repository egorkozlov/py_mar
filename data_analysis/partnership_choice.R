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
directory<-"C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019"
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
ind <-read_dta("C:/Users/Fabio/Dropbox/JMP/empirical analysis/NSFG/88/NSFG88s.dta")


#ind <- ind[!ind$single==1,]
#ind<-ind[ind$state!="California",]
ind$stat<-ind$state
ind$intt<-ind$int1+ind$int3
ind$age2<-ind$age_rel^2
ind$age3<-ind$age_rel^3
#ind<-ind[!is.na(ind$age_rel),]
#ind<-ind[!is.na(ind$birth),]
#Tobit model with individual thresholds: total sample
#for standard errors
G <- length(unique(ind$stat))
c <- G/(G - 1)



lm1 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
           data=ind,#weight=wgt_n,
           na.action=na.exclude)

 
lm2 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
           data=ind[ind$keep==1,],#weight=wgt_n,
           na.action=na.exclude)

lm3 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
           index=c('stat'),
           data=ind[ind$nsfh==1,],#weight=wgt_n,
           na.action=na.exclude)

lm4 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
           index=c('stat'), 
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
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\wmar.tex',sep=""))


#Same but robustess with linear trends
lm1 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind)


lm2 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind[ind$keep==1,])

lm3 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind[ind$nsfh==1,])

lm4 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind[ind$nsfh==0,])


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
              dep.var.labels=c("","","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
              covariate.labels=c("Unilateral Divorce"),
              keep.stat=c("n","rsq"),no.space=TRUE,
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
              keep = c('unid'),
              font.size='footnotesize',
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Linear trend by State", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\wmarlin.tex',sep=""))




##################################################################
##################################################################
####  MARRIAGE versus COHABITATION-compare regimes
##################################################################
##################################################################

ind$intt<-0
ind$intt[ind$tit==0 & ind$unid==1]<-1


#Tobit model with individual thresholds: total sample
#for standard errors

#Summary Statistics
ind1=ind[,c('unid','age_rel','mar','coll','female','birth','nsfh')]

# Rename a column in R
colnames(ind1)[colnames(ind1)=="unid"] <- "Unilateral Divorce Dummy"
colnames(ind1)[colnames(ind1)=="age_rel"] <- "Age Relationship Starts"
colnames(ind1)[colnames(ind1)=="mar"] <- "Married"
colnames(ind1)[colnames(ind1)=="coll"] <- "College"
colnames(ind1)[colnames(ind1)=="female"] <- "Female"
colnames(ind1)[colnames(ind1)=="birth"] <- "Birth year"
colnames(ind1)[colnames(ind1)=="nsfh"] <- "NSFH Dummy"

stat<-stargazer(as.data.frame(ind1),title="Descriptive Statistics: Sample of First Relationships",
                summary.stat=c("N","Mean","Median","sd"),
                float=FALSE,
                out=paste(directory,'\\sum_sing.tex',sep=""))


G <- length(unique(ind$stat))
c <- G/(G - 1)



lm1 <- lm(mar ~ intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind,#weight=wgt_n,
          na.action=na.exclude)


lm2 <- lm(mar ~  intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$keep==1,],#weight=wgt_n,
          na.action=na.exclude)

lm3 <- lm(mar ~  intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$nsfh==1,],#weight=wgt_n,
          na.action=na.exclude)

lm4 <- lm(mar ~  intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$nsfh==0,],#weight=wgt_n,
          na.action=na.exclude)


#Get if coefficients are the same


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


#test_coef_equality(lm1,var1.name = "int2", var2.name = "intt",v=v1)

Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=FALSE,
              se = list(coeftest(lm1,v1)[, 2],
                        coeftest(lm2,v2)[, 2],
                        coeftest(lm3,v3)[, 2],
                        coeftest(lm4,v4)[, 2]),
              omit.table.layout = "n",
              dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
              covariate.labels=c("UnDiv*NoTit","UnDiv*Tit","Tit"),
              keep.stat=c("n","rsq"),no.space=TRUE,
              dep.var.labels=c("","","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              keep = c('intt','int2','tit'),
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),  
              font.size='footnotesize',
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")), 
              out=paste(directory,'\\wmarc.tex',sep=""))


#Same but robustess with linear trends
lm1 <- lm(mar ~ intt+int2+tit+factor(date_rel)+date_rel*stat+age_rel+age2+age3+female,
           data=ind,#weight=wgt_n,
           na.action=na.exclude)


lm2 <- lm(mar ~ intt+int2+tit+factor(date_rel)+date_rel*stat+age_rel+age2+age3+female,
           data=ind[ind$keep==1,],#weight=wgt_n,
           na.action=na.exclude)


lm3 <-lm(mar ~ intt+int2+tit+factor(date_rel)+date_rel*stat+age_rel+age2+age3+female,
         data=ind[ind$nsfh==1,],#weight=wgt_n,
         na.action=na.exclude)

lm4 <- lm(mar ~ intt+int2+tit+factor(date_rel)+date_rel*stat+age_rel+age2+age3+female,
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
              model.names = FALSE,
              dep.var.labels=c("","","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              covariate.labels=c("UnDiv*NoTit","UnDiv*Tit","Tit"),
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
              keep.stat=c("n","rsq"),no.space=TRUE,
              keep = c('intt','int2','tit'),
              font.size='footnotesize',
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Linear trend by State", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\wmarlinc.tex',sep=""))


# ##################################################################
# ##################################################################
# ####  MARRIAGE versus COHABITATION - Logistic Regression
# ##################################################################
# ##################################################################
# 
# #Tobit model with individual thresholds: total sample
# #for standard errors
# G <- length(unique(ind$stat))
# c <- G/(G - 1)
# 
# 
# 
# lm1 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',
#            data=ind,#weight=wgt_n,
#            na.action=na.exclude)
# 
# lm2 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',#weight=wgt_n,
#            data=ind[ind$keep==1,],
#            na.action=na.exclude)
# 
# lm3 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',#weight=wgt_n,
#            data=ind[ind$nsfh==1,],
#            na.action=na.exclude)
# 
# lm4 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',#weight=wgt_n,
#            data=ind[ind$nsfh==0,],
#            na.action=na.exclude)
# 
# 
# #MArginal Effects
# lm1ME<-summary(margins(lm1, variables = "unid"))
# lm2ME<-summary(margins(lm2, variables = "unid"))
# lm3ME<-summary(margins(lm3, variables = "unid"))
# lm4ME<-summary(margins(lm4, variables = "unid"))
# 
# m1<-toString(round(lm1ME$AME[[1]],digits=3))
# m2<-toString(round(lm2ME$AME[[1]],digits=3))
# m3<-toString(round(lm1ME$AME[[1]],digits=3))
# m4<-toString(round(lm2ME$AME[[1]],digits=3))
# 
# #put the stuff in latex
# Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=FALSE,
#               se = list(coeftest(lm1,c*vcovHC(lm1, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm2,c*vcovHC(lm2, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm3,c*vcovHC(lm3, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm4,c*vcovHC(lm4, type = "HC0", cluster = "stat", adjust = T))[, 2]),
#               omit.table.layout = "n",
#               dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
#               covariate.labels=c("Unilateral Divorce"),
#               keep.stat=c("n","rsq"),no.space=TRUE,
#               keep = c('unid'),
#               font.size='footnotesize',
#               column.labels = c("Full Sample", "Resident","NSFH","NSFG"),
#               #column.sep.width = "-15pt",
#               float=FALSE,
#               #float = FALSE,
#               notes = "Standard Errors clustered at the State level.",
#               add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
#                                c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
#                                c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
#                                c("Demographic Controls", "Yes", "Yes", "Yes", "Yes"),
#                                c("\\textbf{Average Marginal Effects}", m1, m2,m3,m4)),
#               out=paste(directory,'\\wmar_logit.tex',sep=""))
# 
# ###########################################################################
# ###########################################################################
# ####  MARRIAGE versus COHABITATION - Logistic Regression with Linear Trends
# ##########################################################################
# ###########################################################################
# G <- length(unique(ind$stat))
# c <- G/(G - 1)
# 
# 
# 
# lm1 <- glm(mar ~ unid+factor(date_rel)+date_rel*stat+age_rel+age2+age3+
#              female+coll+factor(stat)+date_rel*stat,
#            family='binomial',weights=wgt_n,
#            data=ind,
#            na.action=na.exclude)
# 
# lm2 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+
#                  female+coll+factor(stat)+date_rel*stat,
#            family='binomial',weights=wgt_n,
#            data=ind[ind$keep==1,],
#            na.action=na.exclude)
# 
# lm3 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+
#                   female+coll+factor(stat)+date_rel*stat,
#            family='binomial',weights=wgt_n,
#            data=ind[ind$nsfh==1,],
#            na.action=na.exclude)
# 
# lm4 <- glm(mar ~ unid+factor(date_rel)+age_rel+age2+age3+
#                    female+coll+factor(stat)+date_rel*stat,
#            family='binomial',weights=wgt_n,
#            data=ind[ind$nsfh==0,],
#            na.action=na.exclude)
# 
# 
# #MArginal Effects
# lm1ME<-summary(margins(lm1, variables = "unid"))
# lm2ME<-summary(margins(lm2, variables = "unid"))
# lm3ME<-summary(margins(lm3, variables = "unid"))
# lm4ME<-summary(margins(lm4, variables = "unid"))
# 
# m1<-toString(round(lm1ME$AME[[1]],digits=3))
# m2<-toString(round(lm2ME$AME[[1]],digits=3))
# m3<-toString(round(lm1ME$AME[[1]],digits=3))
# m4<-toString(round(lm2ME$AME[[1]],digits=3))
# 
# #put the stuff in latex
# Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=TRUE,
#               se = list(coeftest(lm1,c*vcovHC(lm1, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm2,c*vcovHC(lm2, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm3,c*vcovHC(lm3, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm4,c*vcovHC(lm4, type = "HC0", cluster = "stat", adjust = T))[, 2]),
#               dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
#               covariate.labels=c("Unilateral Divorce"),
#               keep.stat=c("n"),no.space=TRUE,
#               keep = c('unid'),
#               font.size='footnotesize',
#               column.labels = c("Full Sample", "Resident","NSFH","NSFG"),
#               column.sep.width = "-15pt",
#               float = FALSE,
#               notes = "Standard Errors clustered at the State level.",
#               add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
#                                c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
#                                c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
#                                c("Linear trend by State", "Yes", "Yes", "Yes", "Yes"), 
#                                c("Demographic Controls", "Yes", "Yes", "Yes", "Yes"),
#                                c("\\textbf{Average Marginal Effects}", m1, m2,m3,m4)),
#               out='C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019\\wmar_logit_trend.tex')

##############################################################################
##############################################################################
####  MARRIAGE versus COHABITATION - Logistic Regression-HEt by divorce regime
#############################################################################
#############################################################################
# 
# 
# G <- length(unique(ind$stat))
# c <- G/(G - 1)
# 
# 
# 
# lm1 <- glm(mar ~ intt+int2+tit+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',
#            data=ind,
#            na.action=na.exclude)
# 
# lm2 <- glm(mar ~ intt+int2+tit+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',
#            data=ind[ind$keep==1,],
#            na.action=na.exclude)
# 
# lm3 <- glm(mar ~ intt+int2+tit+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',
#            data=ind[ind$nsfh==1,],
#            na.action=na.exclude)
# 
# lm4 <- glm(mar ~ intt+int2+tit+factor(date_rel)+age_rel+age2+age3+female+coll+factor(stat),
#            family='binomial',
#            data=ind[ind$nsfh==0,],
#            na.action=na.exclude)
# 
# 
# #MArginal Effects
# lm1ME<-summary(margins(lm1, variables = c("intt","int2")))
# lm2ME<-summary(margins(lm2, variables = c("intt","int2")))
# lm3ME<-summary(margins(lm3, variables = c("intt","int2")))
# lm4ME<-summary(margins(lm4, variables = c("intt","int2")))
# 
# m1<-toString(round(lm1ME$AME[[1]],digits=3))
# m2<-toString(round(lm2ME$AME[[1]],digits=3))
# m3<-toString(round(lm1ME$AME[[1]],digits=3))
# m4<-toString(round(lm2ME$AME[[1]],digits=3))
# 
# m1b<-toString(round(lm1ME$AME[[2]],digits=3))
# m2b<-toString(round(lm2ME$AME[[2]],digits=3))
# m3b<-toString(round(lm1ME$AME[[2]],digits=3))
# m4b<-toString(round(lm2ME$AME[[2]],digits=3))
# 
# #put the stuff in latex
# Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=TRUE,
#               se = list(coeftest(lm1,c*vcovHC(lm1, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm2,c*vcovHC(lm2, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm3,c*vcovHC(lm3, type = "HC0", cluster = "stat", adjust = T))[, 2],
#                         coeftest(lm4,c*vcovHC(lm4, type = "HC0", cluster = "stat", adjust = T))[, 2]),
#               dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
#               covariate.labels=c("UnDiv*NoTit","UnDiv*Tit","Tit"),
#               keep.stat=c("n"),no.space=TRUE,
#               keep = c('intt','int2','tit'),
#               font.size='footnotesize',
#               column.labels = c("Full Sample", "Resident","NSFH","NSFG"),
#               column.sep.width = "-15pt",
#               float = FALSE,
#               notes = "Standard Errors clustered at the State level.",
#               add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
#                                c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
#                                c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
#                                c("Demographic Controls", "Yes", "Yes", "Yes", "Yes"),
#                                c("\\textbf{AME UnDiv*NoTit}", m1b, m2b,m3b,m4b),
#                                 c("\\textbf{AME UnDiv*Tit}", m1, m2,m3,m4)),
#               out='C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019\\wmar_logit_het.tex')
# 

###########################################
##########################################
#California eliminated
#########################################
##########################################à
ind2<-ind[ind$nsfh==1,]
ind<-ind[ind$state!="California",]

lm1 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind,#weight=wgt_n,
          na.action=na.exclude)


lm2 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$keep==1,],#weight=wgt_n,
          na.action=na.exclude)

lm3 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          index=c('stat'),
          data=ind[ind$nsfh==1,],#weight=wgt_n,
          na.action=na.exclude)

lm4 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
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
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\wmar_c.tex',sep=""))


#Same but robustess with linear trends
lm1 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind)


lm2 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind[ind$keep==1,])

lm3 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind[ind$nsfh==1,])

lm4 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+
            female+coll+stat*date_rel,#weight=wgt_n,
          data=ind[ind$nsfh==0,])


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
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
              keep = c('unid'),
              font.size='footnotesize',
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Linear trend by State", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\wmarlin_c.tex',sep=""))




##################################################################
##################################################################
####  MARRIAGE versus COHABITATION-compare regimes
##################################################################
##################################################################

ind$intt<-0
ind$intt[ind$tit==0 & ind$unid==1]<-1

G <- length(unique(ind$stat))
c <- G/(G - 1)



lm1 <- lm(mar ~ intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind,#weight=wgt_n,
          na.action=na.exclude)


lm2 <- lm(mar ~  intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$keep==1,],#weight=wgt_n,
          na.action=na.exclude)

lm3 <- lm(mar ~  intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$nsfh==1,],#weight=wgt_n,
          na.action=na.exclude)

lm4 <- lm(mar ~  intt+int2+tit+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind[ind$nsfh==0,],#weight=wgt_n,
          na.action=na.exclude)


#Get if coefficients are the same


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
              covariate.labels=c("UnDiv*NoTit","UnDiv*Tit","Tit"),
              dep.var.labels=c("","","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              keep.stat=c("n","rsq"),no.space=TRUE,
              keep = c('intt','int2','tit'),
              column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),  
              font.size='footnotesize',
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"), 
                               c("Birth Year dummies", "Yes", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")), 
              out=paste(directory,'\\wmarc_c.tex',sep=""))


##############################
#Het for children
#############################
lm1 <- lm(mar ~ unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind2[ind2$nch>=0,],#weight=wgt_n,
          na.action=na.exclude)


lm2 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind2[ind2$nch==0,],#weight=wgt_n,
          na.action=na.exclude)

lm3 <- lm(mar ~  unid+factor(date_rel)+factor(birth)+female+coll+factor(stat),
          data=ind2[(ind2$nch==0) & (ind2$nointc==0),],#weight=wgt_n,
          na.action=na.exclude)


#put the stuff in latex
v1<-multiwayvcov::cluster.vcov(lm1, 
                               cbind(ind2[ind2$nch>=0,]$stat, ind2[ind2$nch>=0,]$date_rel),
                               use_white = F, 
                               df_correction = F)

v2<-multiwayvcov::cluster.vcov(lm2, 
                               cbind(ind2[ind2$nch==0,]$stat, ind2[ind2$nch==0,]$date_rel),
                               use_white = F, 
                               df_correction = F)

v3<-multiwayvcov::cluster.vcov(lm3, 
                               cbind(ind2[(ind2$nch==0) & (ind2$nointc==0),]$stat, ind2[(ind2$nch==0) & (ind2$nointc==0),]$date_rel),
                               use_white = F, 
                               df_correction = F)


Df<-stargazer(lm1,lm2,lm3, header=TRUE, align=FALSE,
              se = list(coeftest(lm1,v1)[, 2],
                        coeftest(lm2,v2)[, 2],
                        coeftest(lm3,v3)[, 2]),
              omit.table.layout = "n",
              dep.var.caption=c("\\textit{Dependent variable: Married (0/1)}"),
              dep.var.labels=c("","","\\\\[-4.8ex]"),
              multicolumn = FALSE,
              covariate.labels=c("Unilateral Divorce"),
              keep.stat=c("n","rsq"),no.space=TRUE,
              keep = c('unid'),
              font.size='footnotesize',
              column.labels = c("Some children", "Childless", "Childless+No want children"),
              #column.sep.width = "-15pt",
              float = FALSE,
              #notes = "Standard Errors clustered at the State level.",
              add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes"),
                               c("Year started Fixed Effect", "Yes", "Yes", "Yes"),
                               c("Birth Year dummies", "Yes", "Yes", "Yes"),
                               c("Demographic Controls", "Yes", "Yes", "Yes")),
              out=paste(directory,'\\wmar_chil.tex',sep=""))
