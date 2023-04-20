#get the final sample for the analysis


## remove (almost) everything in the working environment.
rm(list = ls())

#for reading .dta file
library(haven)
library(coda)
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

#fixed effects regression
library(plm)
require(lmtest)

#for exporting in Latex
library(stargazer)
require(foreign)
require(nnet)
require(ggplot2)
require(reshape2)


##################################################################
##################################################################
####  SAMPLE SELECTION AS A FIRST THING
##################################################################
##################################################################
directory<-"C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019"
directory<-"C:\\Users\\Fabio\\Dropbox\\cohabitation and divorce laws\\coh_text"
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
####  COHABITATION LENGTH
##################################################################
##################################################################

ind <-read_dta("C:/Users/Fabio/Dropbox/JMP/empirical analysis/NSFG/88/NSFG88.dta")

#Sample selection
select(ind,16,45,1)

## order by submission_time
ind <- ind[order(ind$t, decreasing=TRUE),]

## remove duplicated user_id
ind <- ind[!duplicated(ind$id),]

#Variable for right censoring
ind$right<-ind$t
ind$right[ind$fail2==1]<-ind$itw[ind$fail2==1]-ind$dati[ind$fail2==1]
ind$right[ind$fail2==2]<-ind$itw[ind$fail2==2]-ind$dati[ind$fail2==2]


#Drop if inconsistency
ind<-ind[ind$right>=ind$t,]
ind$age<-ind$agep
ind$age2<-ind$age^2
ind$age3<-ind$age^3
ind$Married<-0
ind$Married[ind$fail2==1]<-1
ind$Separated<-0
ind$Separated[ind$fail2==2]<-1
ind$Married1<-1
ind$Married1[ind$fail2==1]<-2
ind$Separated1<-1
ind$Separated1[ind$fail2==2]<-2
#Summary Statistics
ind1=ind[,c('unid','age','beg','coll','female','t','birth','nsfh','censored','Married','Separated')]

# Rename a column in R
colnames(ind1)[colnames(ind1)=="unid"] <- "Unilateral Divorce"
colnames(ind1)[colnames(ind1)=="age"] <- "Age Cohabitation Starts"
colnames(ind1)[colnames(ind1)=="beg"] <- "Year Cohabitation Starts"
colnames(ind1)[colnames(ind1)=="coll"] <- "College"
colnames(ind1)[colnames(ind1)=="female"] <- "Female"
colnames(ind1)[colnames(ind1)=="t"] <- "Cohabitation Duration"
colnames(ind1)[colnames(ind1)=="birth"] <- "Year of birth"
colnames(ind1)[colnames(ind1)=="nsfh"] <- "NSFH"
colnames(ind1)[colnames(ind1)=="censored"] <- "Censored"

stat<-stargazer(as.data.frame(ind1),title="Descriptive Statistics: Sample of First Cohabitation",
                summary.stat=c("N","Mean","Median","sd"),
                float=FALSE,
                out=paste(directory,'\\sum_coh.tex',sep=""))

res.cox <- coxph(Surv(t, Married1) ~ unid+factor(stat)+age+age2+age3+factor(beg)+coll+female, data = ind)

res.cox1 <- coxph(Surv(t, Separated1) ~ unid+factor(stat)+age+age2+age3+factor(beg)+coll+female, data = ind)

# #Tobit model with individual thresholds: total sample
# tob1 <- vglm(t ~ unid+factor(stat)+age+age2+age3+factor(beg)+coll+female,
#              tobit(Lower = -Inf, Upper = ind$right,imethod = 2),
#              data = ind)
# 
# 
# 
# #Tobit model with individual thresholds: only always leaved there
# tob2 <- vglm(t ~ unid+factor(stat)+age+age2+age3+factor(beg)+coll+female,
#              tobit(Lower = -Inf, Upper = ind$right[ind$keep==1],imethod = 2),
#              data = ind[ind$keep==1,])
# 
# #Tobit model with individual thresholds: total sample
# tob3 <- vglm(t ~ unid+factor(stat)+age+age2+age3+factor(beg)+coll+female,
#              tobit(Lower = -Inf, Upper = ind$right[ind$nsfh==1],imethod = 2),
#              data = ind[ind$nsfh==1,])
# 
# 
# 
# #Tobit model with individual thresholds: total sample
# tob4 <- vglm(t ~ unid+factor(stat)+age+age2+age3+factor(beg)+coll,
#              tobit(Lower = -Inf, Upper = ind$right[ind$nsfh==0],imethod = 2),
#              data = ind[ind$nsfh==0,])
# 
# #Marginal Effects(NOT on the latent variable...)
# IMR <- Vectorize( function(x) exp( dnorm(x,log=T) - pnorm(x,log.p = T) ) )
# #Tob1
# 
# mewdn <- ind
# mewdn$unid<-0
# mewdy <- ind
# mewdy$unid<-1
# fitn <- predict(tob1, newdata = mewdn, type = "response")
# fity <- predict(tob1, newdata = mewdy, type = "response")
# av1<-mean(fity-fitn)*pnorm(mean(predict(tob1))/exp(coef(tob1)[2]))
# av1<-mean(fity+exp(coef(tob1)[2])*IMR(((-fity+ind$right)/exp(coef(tob1)[2])))-
#             (fitn+exp(coef(tob1)[2])*IMR(((-fitn+ind$right)/exp(coef(tob1)[2])))))
# 
# 
# #Tob2
# 
# mewdn <- ind[ind$keep==1,]
# mewdn$unid<-0
# mewdy <- ind[ind$keep==1,]
# mewdy$unid<-1
# fitn <- predict(tob2, newdata = mewdn, type = "response")
# fity <- predict(tob2, newdata = mewdy, type = "response")
# av2<-mean(fity-fitn)*pnorm(mean(predict(tob2))/exp(coef(tob2)[2]))
# av2<-mean(fity+exp(coef(tob2)[2])*IMR(((-fity+ind$right[ind$keep==1])/exp(coef(tob2)[2])))-
#             (fitn+exp(coef(tob2)[2])*IMR(((-fitn+ind$right[ind$keep==1])/exp(coef(tob2)[2])))))
# 
# 
# #Tob3
# 
# mewdn <- ind[ind$nsfh==1,]
# mewdn$unid<-0
# mewdy <- ind[ind$nsfh==1,]
# mewdy$unid<-1
# fitn <- predict(tob3, newdata = mewdn, type = "response")
# fity <- predict(tob3, newdata = mewdy, type = "response")
# av3<-mean(fity-fitn)*pnorm(mean(predict(tob3))/exp(coef(tob3)[2]))
# av3<-mean(fity+exp(coef(tob3)[2])*IMR(((-fity+ind$right[ind$nsfh==1])/exp(coef(tob3)[2])))-
#             (fitn+exp(coef(tob3)[2])*IMR(((-fitn+ind$right[ind$nsfh==1])/exp(coef(tob3)[2])))))
# 
# 
# #Tob4
# 
# mewdn <- ind[ind$nsfh==0,]
# mewdn$unid<-0
# mewdy <- ind[ind$nsfh==0,]
# mewdy$unid<-1
# fitn <- predict(tob4, newdata = mewdn, type = "response")
# fity <- predict(tob4, newdata = mewdy, type = "response")
# av4<-mean(fity-fitn)*pnorm(mean(predict(tob4))/exp(coef(tob4)[2]))
# av4<-mean(fity+exp(coef(tob4)[2])*IMR(((-fity+ind$right[ind$nsfh==0])/exp(coef(tob4)[2])))-
#             (fitn+exp(coef(tob4)[2])*IMR(((-fitn+ind$right[ind$nsfh==0])/exp(coef(tob4)[2])))))
# 
# #Put Stuff in a table
# #Build the tables for the paper
# sink("temp.Rnw")
# 
# cat(paste(
#   '\\begin{table}[!htbp] \\centering 
# \\scriptsize
# \\begin{tabular}{@{\\extracolsep{5pt}}lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} } 
# \\\\[-1.8ex]\\hline 
# \\hline \\\\[-1.8ex] 
#  & \\multicolumn{4}{c}{\\textit{Dependent variable: Cohabitation Length}} \\\\ 
# \\cline{2-5}
# \\\\[-1.8ex] & \\multicolumn{1}{c}{Full Sample} & \\multicolumn{1}{c}{Resident}& \\multicolumn{1}{c}{NSFH}& \\multicolumn{1}{c}{NSFG} \\\\ 
# \\\\[-1.8ex] & \\multicolumn{1}{c}{(1)} & \\multicolumn{1}{c}{(2)} & \\multicolumn{1}{c}{(3)} & \\multicolumn{1}{c}{(4)}\\\\ 
# \\hline \\\\[-1.8ex] 
#  Unilateral Divorce & ',round(coef(tob1)[3],digits=2),'^{***} & ',round(coef(tob2)[3],digits=2),'^{***} & ',round(coef(tob3)[3],digits=2),'^{***} & ',round(coef(tob4)[3],digits=2),'^{***} \\\\ 
#   & (',round(sqrt(vcov(tob1)[3,3]),digits=2),') & (',round(sqrt(vcov(tob2)[3,3]),digits=2),') & (',round(sqrt(vcov(tob3)[3,3]),digits=2),') & (',round(sqrt(vcov(tob4)[3,3]),digits=2),') \\\\  
#  \\hline \\\\[-1.8ex]
# Marginal Effect & ',round(av1,digits=2),' & ',round(av2,digits=2),' & ',round(av3,digits=2),' & ',round(av4,digits=2),' \\\\ 
#  \\hline \\\\[-1.8ex] 
# State Fixed effects & Yes & Yes & Yes & Yes \\\\ 
# Year Fixed effects & Yes & Yes & Yes & Yes \\\\ 
# Age Polynomial & Yes & Yes & Yes & Yes \\\\ 
# Demographic Controls & Yes & Yes & Yes & Yes \\\\ 
# Observations & \\multicolumn{1}{c}{',nobs(tob1),'} & \\multicolumn{1}{c}{',nobs(tob2),'} & \\multicolumn{1}{c}{',nobs(tob3),'} & \\multicolumn{1}{c}{',nobs(tob4),'} \\\\ 
# \\hline
# Censored(\\%) & \\multicolumn{1}{c}{',round(mean(ind$censored),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$keep==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$nsfh==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$nsfh==0]),digits=2),'} \\\\ 
# \\hline 
# \\hline \\\\[-1.8ex] 
# \\textit{Note:}  & \\multicolumn{4}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\\\ 
# \\end{tabular} 
# \\end{table} '))
# 
# sink()
# Sweave("temp.Rnw")
# 
# 
# #Robustness: OLS with and without time trends
# #Tobit model with individual thresholds: total sample
# G <- length(unique(ind$stat))
# c <- G/(G - 1)
# 
# lm1 <- plm(t ~ unid+factor(beg)+age+age2+age3+coll+female,
#            index=c('stat'), model = 'within',
#            data=ind,
#            na.action=na.exclude,
#            effect='individual')
# 
# lm2 <- plm(t ~ unid+factor(beg)+age+age2+age3+coll+female,
#            index=c('stat'), model = 'within',
#            data=ind[ind$keep==1,],
#            na.action=na.exclude,
#            effect='individual')
# 
# #Drop Censored
# ind<-ind[ind$right>ind$t,]
# 
# lm3 <- plm(t ~ unid+factor(beg)+age+age2+age3+coll+female,
#            index=c('stat'), model = 'within',
#            data=ind[ind$nsfh==1,],
#            na.action=na.exclude,
#            effect='individual')
# 
# lm4 <- plm(t ~ unid+factor(beg)+age+age2+age3+coll,
#            index=c('stat'), model = 'within',
#            data=ind[ind$nsfh==0,],
#            na.action=na.exclude,
#            effect='individual')
# 
# Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=TRUE,
#               se = list(coeftest(lm1,c*vcovHC(lm1, type = "HC1", cluster = "group"))[, 2],
#                         coeftest(lm2,c*vcovHC(lm2, type = "HC1", cluster = "group"))[, 2],
#                         coeftest(lm3,c*vcovHC(lm3, type = "HC1", cluster = "group"))[, 2],
#                         coeftest(lm4,c*vcovHC(lm4, type = "HC1", cluster = "group"))[, 2]),
#               dep.var.labels=c("Marriage Length"),
#               column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
#               float = FALSE,
#               covariate.labels=c("Unilateral Divorce"),
#               keep.stat=c("n"),no.space=TRUE,
#               keep = c('unid'),
#               font.size='footnotesize',
#               column.sep.width = "-15pt",
#               notes = "Standard Errors clustered at the State level.",
#               add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
#                                c("Age Polynomials", "Yes", "Yes", "Yes", "Yes"),
#                                c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
#                                c("Linear trend by State", "No", "No", "No", "No"),
#                                c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
#               out='C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019\\ols.tex')
# 
# 
# #Robustness: OLS with and without time trends
# #Tobit model with individual thresholds: total sample
# G <- length(unique(ind$stat))
# c <- G/(G - 1)
# 
# lm1 <- plm(t ~ unid+factor(beg)+factor(birth)+beg*stat+coll+female,
#            index=c('stat'), model = 'within',
#            data=ind,
#            na.action=na.exclude,
#            effect='individual')
# 
# lm2 <- plm(t ~ unid+factor(beg)+age+age2+age3+beg*stat+coll+female,
#            index=c('stat'), model = 'within',
#            data=ind[ind$keep==1,],
#            na.action=na.exclude,
#            effect='individual')
# 
# #Drop Censored
# ind<-ind[ind$right>ind$t,]
# 
# lm3 <- plm(t ~ unid+factor(beg)+age+age2+age3+beg*stat+coll+female,
#            index=c('stat'), model = 'within',
#            data=ind[ind$nsfh==1,],
#            na.action=na.exclude,
#            effect='individual')
# 
# lm4 <- plm(t ~ unid+factor(beg)+age+age2+age3+beg*stat+coll,
#            index=c('stat'), model = 'within',
#            data=ind[ind$nsfh==0,],
#            na.action=na.exclude,
#            effect='individual')
# 
# Df<-stargazer(lm1,lm2,lm3,lm4, header=TRUE, align=TRUE,
#               se = list(coeftest(lm1,c*vcovHC(lm1, type = "HC1", cluster = "group"))[, 2],
#                         coeftest(lm2,c*vcovHC(lm2, type = "HC1", cluster = "group"))[, 2],
#                         coeftest(lm3,c*vcovHC(lm3, type = "HC1", cluster = "group"))[, 2],
#                         coeftest(lm4,c*vcovHC(lm4, type = "HC1", cluster = "group"))[, 2]),
#               dep.var.labels=c("Marriage Length"),
#               covariate.labels=c("Unilateral Divorce"),
#               column.labels = c("Full Sample", "Resident", "NSFH", "NSFG"),
#               float = FALSE,
#               keep.stat=c("n"),no.space=TRUE,
#               keep = c('unid'),
#               font.size='footnotesize',
#               column.sep.width = "-15pt",
#               label="table:reg",
#               notes = "Standard Errors clustered at the State level.",
#               add.lines = list(c("State Fixed effects", "Yes", "Yes", "Yes", "Yes"),
#                                c("Age Polynomials", "Yes", "Yes", "Yes", "Yes"),
#                                c("Year started Fixed Effect", "Yes", "Yes", "Yes", "Yes"),
#                                c("Linear trend by State", "Yes", "Yes", "Yes", "Yes"),
#                                c("Demographic Controls", "Yes", "Yes", "Yes", "Yes")),
#               out='C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019\\olslin.tex')
# 
# 
# #Heterogeneity By divorce regime
# #Tobit model with individual thresholds: total sample
# ind$intt<-0
# ind$intt[ind$tit==0 & ind$unid==1]<-1
# 
# tob1 <- vglm(t ~ intt+int2+tit+factor(stat)+age+age2+age3+factor(beg)+coll+female,
#              tobit(Lower = -Inf, Upper = ind$right,imethod = 2),
#              data = ind)
# 
# 
# 
# #Tobit model with individual thresholds: only always leaved there
# tob2 <- vglm(t ~ intt+int2+tit+factor(stat)+age+age2+age3+factor(beg)+coll+female,
#              tobit(Lower = -Inf, Upper = ind$right[ind$keep==1],imethod = 2),
#              data = ind[ind$keep==1,])
# 
# #Tobit model with individual thresholds: total sample
# tob3 <- vglm(t ~intt+int2+tit+factor(stat)+age+age2+age3+factor(beg)+coll+female,
#              tobit(Lower = -Inf, Upper = ind$right[ind$nsfh==1],imethod = 2),
#              data = ind[ind$nsfh==1,])
# 
# 
# 
# #Tobit model with individual thresholds: total sample
# tob4 <- vglm(t ~ intt+int2+tit+factor(stat)+age+age2+age3+factor(beg)+coll,
#              tobit(Lower = -Inf, Upper = ind$right[ind$nsfh==0],imethod = 2),
#              data = ind[ind$nsfh==0,])
# 
# #Marginal Effects(NOT on the latent variable...)
# IMR <- Vectorize( function(x) exp( dnorm(x,log=T) - pnorm(x,log.p = T) ) )
# #Tob1
# 
# mewdn <- ind
# mewdn$intt<-0
# mewdy <- ind
# mewdy$intt<-1
# fitn <- predict(tob1, newdata = mewdn, type = "response")
# fity <- predict(tob1, newdata = mewdy, type = "response")
# av1<-mean(fity-fitn)*pnorm(mean(predict(tob1))/exp(coef(tob1)[3]))
# av1<-mean(fity+exp(coef(tob1)[3])*IMR(((-fity+ind$right)/exp(coef(tob1)[3])))-
#             (fitn+exp(coef(tob1)[3])*IMR(((-fitn+ind$right)/exp(coef(tob1)[3])))))
# 
# 
# #Tob2
# 
# mewdn <- ind[ind$keep==1,]
# mewdn$intt<-0
# mewdy <- ind[ind$keep==1,]
# mewdy$intt<-1
# fitn <- predict(tob2, newdata = mewdn, type = "response")
# fity <- predict(tob2, newdata = mewdy, type = "response")
# av2<-mean(fity-fitn)*pnorm(mean(predict(tob2))/exp(coef(tob2)[3]))
# av2<-mean(fity+exp(coef(tob2)[3])*IMR(((-fity+ind$right[ind$keep==1])/exp(coef(tob2)[3])))-
#             (fitn+exp(coef(tob2)[3])*IMR(((-fitn+ind$right[ind$keep==1])/exp(coef(tob2)[3])))))
# 
# 
# #Tob3
# 
# mewdn <- ind[ind$nsfh==1,]
# mewdn$intt<-0
# mewdy <- ind[ind$nsfh==1,]
# mewdy$intt<-1
# fitn <- predict(tob3, newdata = mewdn, type = "response")
# fity <- predict(tob3, newdata = mewdy, type = "response")
# av3<-mean(fity-fitn)*pnorm(mean(predict(tob3))/exp(coef(tob3)[2]))
# av3<-mean(fity+exp(coef(tob3)[3])*IMR(((-fity+ind$right[ind$nsfh==1])/exp(coef(tob3)[3])))-
#             (fitn+exp(coef(tob3)[3])*IMR(((-fitn+ind$right[ind$nsfh==1])/exp(coef(tob3)[3])))))
# 
# 
# #Tob4
# 
# mewdn <- ind[ind$nsfh==0,]
# mewdn$intt<-0
# mewdy <- ind[ind$nsfh==0,]
# mewdy$intt<-1
# fitn <- predict(tob4, newdata = mewdn, type = "response")
# fity <- predict(tob4, newdata = mewdy, type = "response")
# av4<-mean(fity-fitn)*pnorm(mean(predict(tob4))/exp(coef(tob4)[3]))
# av4<-mean(fity+exp(coef(tob4)[3])*IMR(((-fity+ind$right[ind$nsfh==0])/exp(coef(tob4)[3])))-
#             (fitn+exp(coef(tob4)[3])*IMR(((-fitn+ind$right[ind$nsfh==0])/exp(coef(tob4)[3])))))
# 
# 
# #Tob1
# 
# mewdn <- ind
# mewdn$int2<-0
# mewdy <- ind
# mewdy$int2<-1
# fitn <- predict(tob1, newdata = mewdn, type = "response")
# fity <- predict(tob1, newdata = mewdy, type = "response")
# av1a<-mean(fity-fitn)*pnorm(mean(predict(tob1))/exp(coef(tob1)[4]))
# av1a<-mean(fity+exp(coef(tob1)[4])*IMR(((-fity+ind$right)/exp(coef(tob1)[4])))-
#              (fitn+exp(coef(tob1)[4])*IMR(((-fitn+ind$right)/exp(coef(tob1)[4])))))
# 
# 
# #Tob2
# 
# mewdn <- ind[ind$keep==1,]
# mewdn$int2<-0
# mewdy <- ind[ind$keep==1,]
# mewdy$int2<-1
# fitn <- predict(tob2, newdata = mewdn, type = "response")
# fity <- predict(tob2, newdata = mewdy, type = "response")
# av2a<-mean(fity-fitn)*pnorm(mean(predict(tob2))/exp(coef(tob2)[4]))
# av2a<-mean(fity+exp(coef(tob2)[4])*IMR(((-fity+ind$right[ind$keep==1])/exp(coef(tob2)[4])))-
#              (fitn+exp(coef(tob2)[4])*IMR(((-fitn+ind$right[ind$keep==1])/exp(coef(tob2)[4])))))
# 
# 
# #Tob3
# 
# mewdn <- ind[ind$nsfh==1,]
# mewdn$int2<-0
# mewdy <- ind[ind$nsfh==1,]
# mewdy$int2<-1
# fitn <- predict(tob3, newdata = mewdn, type = "response")
# fity <- predict(tob3, newdata = mewdy, type = "response")
# av3a<-mean(fity-fitn)*pnorm(mean(predict(tob3))/exp(coef(tob3)[2]))
# av3a<-mean(fity+exp(coef(tob3)[4])*IMR(((-fity+ind$right[ind$nsfh==1])/exp(coef(tob3)[4])))-
#              (fitn+exp(coef(tob3)[4])*IMR(((-fitn+ind$right[ind$nsfh==1])/exp(coef(tob3)[4])))))
# 
# 
# #Tob4
# 
# mewdn <- ind[ind$nsfh==0,]
# mewdn$int2<-0
# mewdy <- ind[ind$nsfh==0,]
# mewdy$int2<-1
# fitn <- predict(tob4, newdata = mewdn, type = "response")
# fity <- predict(tob4, newdata = mewdy, type = "response")
# av4a<-mean(fity-fitn)*pnorm(mean(predict(tob4))/exp(coef(tob4)[4]))
# av4a<-mean(fity+exp(coef(tob4)[4])*IMR(((-fity+ind$right[ind$nsfh==0])/exp(coef(tob4)[4])))-
#              (fitn+exp(coef(tob4)[4])*IMR(((-fitn+ind$right[ind$nsfh==0])/exp(coef(tob4)[4])))))
# 
# 
# 
# #Put Stuff in a table
# #Build the tables for the paper
# sink("temphet.Rnw")
# 
# cat(paste(
#   '\\begin{table}[!htbp] \\centering 
#   \\label{table:reg} 
# \\scriptsize
# \\begin{tabular}{@{\\extracolsep{5pt}}lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} } 
# \\\\[-1.8ex]\\hline 
# \\hline \\\\[-1.8ex] 
#  & \\multicolumn{4}{c}{\\textit{Dependent variable: Cohabitation Length}} \\\\ 
# \\cline{2-5}
# \\\\[-1.8ex] & \\multicolumn{1}{c}{Full Sample} & \\multicolumn{1}{c}{Resident}& \\multicolumn{1}{c}{NSFH}& \\multicolumn{1}{c}{NSFG} \\\\ 
# \\\\[-1.8ex] & \\multicolumn{1}{c}{(1)} & \\multicolumn{1}{c}{(2)} & \\multicolumn{1}{c}{(3)} & \\multicolumn{1}{c}{(4)}\\\\ 
# \\hline \\\\[-1.8ex] 
#  UniDiv*NoTit & ',round(coef(tob1)[3],digits=2),'^{***} & ',round(coef(tob2)[3],digits=2),'^{***} & ',round(coef(tob3)[3],digits=2),'^{***} & ',round(coef(tob4)[3],digits=2),'^{***} \\\\ 
#   & (',round(sqrt(vcov(tob1)[3,3]),digits=2),') & (',round(sqrt(vcov(tob2)[3,3]),digits=2),') & (',round(sqrt(vcov(tob3)[3,3]),digits=2),') & (',round(sqrt(vcov(tob4)[3,3]),digits=2),') \\\\  
#  \\hline \\\\[-1.8ex]
#  UniDiv*Tit & ',round(coef(tob1)[4],digits=2),'^{***} & ',round(coef(tob2)[4],digits=2),'^{***} & ',round(coef(tob3)[4],digits=2),'^{***} & ',round(coef(tob4)[4],digits=2),'^{***} \\\\ 
#   & (',round(sqrt(vcov(tob1)[4,4]),digits=2),') & (',round(sqrt(vcov(tob2)[4,4]),digits=2),') & (',round(sqrt(vcov(tob3)[4,4]),digits=2),') & (',round(sqrt(vcov(tob4)[4,4]),digits=2),') \\\\  
#  \\hline \\\\[-1.8ex]
#   Tit & ',round(coef(tob1)[5],digits=2),'^{***} & ',round(coef(tob2)[5],digits=2),'^{***} & ',round(coef(tob3)[5],digits=2),'^{***} & ',round(coef(tob4)[5],digits=2),'^{***} \\\\ 
#   & (',round(sqrt(vcov(tob1)[5,5]),digits=2),') & (',round(sqrt(vcov(tob2)[5,5]),digits=2),') & (',round(sqrt(vcov(tob3)[5,5]),digits=2),') & (',round(sqrt(vcov(tob4)[5,5]),digits=2),') \\\\  
#  \\hline \\\\[-1.8ex]
# State Fixed effects & Yes & Yes & Yes & Yes \\\\ 
# Year Fixed effects & Yes & Yes & Yes & Yes \\\\ 
# Age Polynomial & Yes & Yes & Yes & Yes \\\\ 
# Demographic Controls & Yes & Yes & Yes & Yes \\\\ 
# Observations & \\multicolumn{1}{c}{',nobs(tob1),'} & \\multicolumn{1}{c}{',nobs(tob2),'} & \\multicolumn{1}{c}{',nobs(tob3),'} & \\multicolumn{1}{c}{',nobs(tob4),'} \\\\ 
# \\hline
# Censored(\\%) & \\multicolumn{1}{c}{',round(mean(ind$censored),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$keep==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$nsfh==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$nsfh==0]),digits=2),'} \\\\ 
# \\hline 
# \\hline \\\\[-1.8ex] 
# \\textit{Note:}  & \\multicolumn{4}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\\\ 
# \\end{tabular} 
# \\end{table} '))
# 
# sink()
# Sweave("temphet.Rnw")
# ##################################################################
# ##################################################################
####  COHABITATION END-MULTINOMIAL LOGIT
##################################################################
##################################################################
#Interpretation:https://stats.idre.ucla.edu/r/dae/multinomial-logistic-regression/

ind <-read_dta("C:/Users/Fabio/Dropbox/JMP/empirical analysis/NSFG/88/NSFG88.dta")

#Sample selection
select(ind,16,45,1)
ind$age<-ind$agep
ind$age2<-ind$age^2
ind$age3<-ind$age^3

#Reshape Data
tm<-mlogit.data(ind,choice="fail2",shape="wide",id.var = "nid")
# 
# #MLOGIT model with individual thresholds: total sample
# ml1 <- mnlogit(fail2 ~ 0 |  unid+factor(stat)+factor(tt)+factor(beg)+age+age2,
#                data=tm,reflevel="0")
# 
# #mLOGIT model with individual thresholds: only always leaved there
# ml2 <- mnlogit(fail2 ~ 0 | unid+coll+factor(stat)+factor(tt)+age+age2+age3+factor(beg)+female+coll,
#                data =tm[tm$keep==1,],reflevel="0")
# 
# #MLOGIT model with individual thresholds: NSFG
# ml3 <- mnlogit(fail2 ~ 0 | unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg)+female+coll,
#                data =tm[tm$nsfh==1,],reflevel="0")
# 
# 
# #MLOGIT model with individual thresholds: NSFG
# ml4 <- mnlogit(fail2 ~ 0 | unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg)+coll,
#                data =tm[tm$nsfh==0,],reflevel="0")
# 
# #Get the effects
# p1a<-coef(ml1)[grep("unid:1$", names(coef(ml1)))]
# s1a<-sqrt(vcov(ml1)[grep("unid:1$", names(coef(ml1))),grep("unid:1$", names(coef(ml1)))])
# p1b<-coef(ml1)[grep("unid:2$", names(coef(ml1)))]
# s1b<-sqrt(vcov(ml1)[grep("unid:2$", names(coef(ml1))),grep("unid:2$", names(coef(ml1)))])
# 
# p2a<-coef(ml2)[grep("unid:1$", names(coef(ml2)))]
# s2a<-sqrt(vcov(ml2)[grep("unid:1$", names(coef(ml2))),grep("unid:1$", names(coef(ml2)))])
# p2b<-coef(ml2)[grep("unid:2$", names(coef(ml2)))]
# s2b<-sqrt(vcov(ml2)[grep("unid:2$", names(coef(ml2))),grep("unid:2$", names(coef(ml2)))])
# 
# p3a<-coef(ml3)[grep("unid:1$", names(coef(ml3)))]
# s3a<-sqrt(vcov(ml3)[grep("unid:1$", names(coef(ml3))),grep("unid:1$", names(coef(ml3)))])
# p3b<-coef(ml3)[grep("unid:2$", names(coef(ml3)))]
# s3b<-sqrt(vcov(ml3)[grep("unid:2$", names(coef(ml3))),grep("unid:2$", names(coef(ml3)))])
# 
# p4a<-coef(ml4)[grep("unid:1$", names(coef(ml4)))]
# s4a<-sqrt(vcov(ml4)[grep("unid:1$", names(coef(ml4))),grep("unid:1$", names(coef(ml4)))])
# p4b<-coef(ml4)[grep("unid:2$", names(coef(ml4)))]
# s4b<-sqrt(vcov(ml4)[grep("unid:2$", names(coef(ml4))),grep("unid:2$", names(coef(ml4)))])
# 
# 
# #Put in a table
# sink("mlc.Rnw")
# 
# cat(paste(
#   '\\begin{table}[!htbp] \\centering 
#   \\caption{} 
#   \\label{table:reg} 
# \\footnotesize 
# \\begin{tabular}{@{\\extracolsep{5pt}}lD{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} D{.}{.}{-3} } 
# \\\\[-1.8ex]\\hline 
# \\hline 
# \\\\[-1.8ex] & \\multicolumn{1}{c}{Full Sample} & \\multicolumn{1}{c}{Resident}& \\multicolumn{1}{c}{NSFH}& \\multicolumn{1}{c}{NSFG} \\\\ 
# \\\\[-1.8ex] & \\multicolumn{1}{c}{(1)} & \\multicolumn{1}{c}{(2)} & \\multicolumn{1}{c}{(3)} & \\multicolumn{1}{c}{(4)}\\\\ 
# \\hline \\\\[-1.8ex] 
# \\\\[-2.2ex] & \\multicolumn{4}{c}{\\textit{Dependent variable: $\\ \\log(\\text{Pr Married}/\\text{Pr Cohabit})$}} \\\\  
#  \\hline \\\\[-1.8ex]
#  Unilateral Divorce & ',round(p1a,digits=2),' & ',round(p2a,digits=2),' & ',round(p3a,digits=2),' & ',round(p4a,digits=2),' \\\\ 
#   & (',round(s1a,digits=2),') & (',round(s2a,digits=2),') & (',round(s3a,digits=2),') & (',round(s4a,digits=2),') \\\\  
#  \\hline \\\\[-1.8ex]
#  \\\\[-2.2ex] & \\multicolumn{4}{c}{\\textit{Dependent variable: $\\ \\log(\\text{Pr Separate}/\\text{Pr Cohabit})$}} \\\\  
#  \\hline \\\\[-1.8ex]
#  Unilateral Divorce & ',round(p1b,digits=2),' & ',round(p2b,digits=2),' & ',round(p3b,digits=2),' & ',round(p4b,digits=2),' \\\\ 
#   & (',round(s1b,digits=2),') & (',round(s2b,digits=2),') & (',round(s3b,digits=2),') & (',round(s4b,digits=2),') \\\\  
#  \\hline \\\\[-1.8ex]
# State Fixed effects & Yes & Yes & Yes & Yes \\\\ 
# Year Fixed effects & Yes & Yes & Yes & Yes \\\\ 
# Age Polynomial & Yes & Yes & Yes & Yes \\\\
# Picewise Duration & Yes & Yes & Yes & Yes \\\\ 
# Demographic Controls & Yes & Yes & Yes & Yes \\\\ 
# \\hline
# Observations & \\multicolumn{1}{c}{',round(length(ind$id)/3,digits=2),'} & \\multicolumn{1}{c}{',round(length(ind$id[ind$keep==1])/3,digits=2),'} & \\multicolumn{1}{c}{',round(length(ind$id[ind$nsfh==1])/3,digits=2),'} & \\multicolumn{1}{c}{',round(length(ind$id[ind$nsfh==0])/3,digits=2),'} \\\\ 
# \\hline
# Individuals & \\multicolumn{1}{c}{',nobs(tob1),'} & \\multicolumn{1}{c}{',nobs(tob2),'} & \\multicolumn{1}{c}{',nobs(tob3),'} & \\multicolumn{1}{c}{',nobs(tob4),'} \\\\ 
# \\hline
# Censored Individuals(\\%) & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$sel==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$keep==1 & ind$sel==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$nsfh==1 & ind$sel==1]),digits=2),'} & \\multicolumn{1}{c}{',round(mean(ind$censored[ind$nsfh==0 & ind$sel==1]),digits=2),'} \\\\ 
# \\hline 
# \\hline \\\\[-1.8ex] 
# \\textit{Note:}  & \\multicolumn{4}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\\\ 
# \\end{tabular} 
# \\end{table} '))
# 
# sink()
# Sweave("mlc.Rnw")
# 
# #############$$$$$$$$$$$$$$$$$
##HMultinomial
#############$$$$$$$$$$$$$$$$
#Mprobit model with individual thresholds: total sample
# ml1 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+factor(beg)+age+age2+age3,
#            data =ind, n.draws = 50, burnin = 10,thin = 3, verbose = TRUE)
# 
# #mprobit model with individual thresholds: only always leaved there
# ml2 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg),
#            data =ind[ind$keep==1,], n.draws = 50, burnin = 10,thin = 3, verbose = TRUE)
# 
# #Mprobit model with individual thresholds: NSFG
# ml3 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg),
#            data =ind[ind$nsfh==1,], n.draws = 50, burnin = 10,thin = 3, verbose = TRUE)
# 
# 
# #Mprobit model with individual thresholds: NSFG
# ml4 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg),
#            data =ind[ind$nsfh==0,], n.draws = 500, burnin = 10,thin = 3, verbose = TRUE)

ml1 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+factor(beg)+age+age2+age3,
           data =ind, n.draws = 50,  verbose = TRUE)

#mprobit model with individual thresholds: only always leaved there
ml2 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg),
           data =ind[ind$keep==1,], n.draws = 50,  verbose = TRUE)

#Mprobit model with individual thresholds: NSFG
ml3 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg),
           data =ind[ind$nsfh==1,], n.draws = 50,  verbose = TRUE)


#Mprobit model with individual thresholds: NSFG
ml4 <- mnp(fail2 ~  unid+factor(stat)+factor(tt)+age+age2+age3+factor(beg),
           data =ind[ind$nsfh==0,], n.draws = 50,  verbose = TRUE)







#Get the effects
p1a<-summary(ml1)$coef[3,1]
s1a<-summary(ml1)$coef[3,2]
p1b<-summary(ml1)$coef[4,1]
s1b<-summary(ml1)$coef[4,2]

p2a<-summary(ml2)$coef[3,1]
s2a<-summary(ml2)$coef[3,2]
p2b<-summary(ml2)$coef[4,1]
s2b<-summary(ml2)$coef[4,2]

p3a<-summary(ml3)$coef[3,1]
s3a<-summary(ml3)$coef[3,2]
p3b<-summary(ml3)$coef[4,1]
s3b<-summary(ml3)$coef[4,2]

p4a<-summary(ml4)$coef[3,1]
s4a<-summary(ml4)$coef[3,2]
p4b<-summary(ml4)$coef[4,1]
s4b<-summary(ml4)$coef[4,2]

#Make a comparable coefffficient to mlogit
mewdn <- ind
mewdn$unid<-0
mewdy <- ind
mewdy$unid<-1
fitnc1<-mean(predict(ml1, newdata = mewdn, type = "prob")$p[,1])
fityc1 <- mean(predict(ml1, newdata = mewdy, type = "prob")$p[,1])
fitnm1<-mean(predict(ml1, newdata = mewdn, type = "prob")$p[,2])
fitym1 <- mean(predict(ml1, newdata = mewdy, type = "prob")$p[,2])
fitns1<-mean(predict(ml1, newdata = mewdn, type = "prob")$p[,3])
fitys1<- mean(predict(ml1, newdata = mewdy, type = "prob")$p[,3])
efm1<-exp(log(fitym1/fityc1)-log(fitnm1/fitnc1))
efs1<-exp(log(fitys1/fityc1)-log(fitns1/fitnc1))

mewdn <- ind[ind$keep==1,]
mewdn$unid<-0
mewdy <- ind[ind$keep==1,]
mewdy$unid<-1
fitnc2<-mean(predict(ml2, newdata = mewdn, type = "prob")$p[,1])
fityc2 <- mean(predict(ml2, newdata = mewdy, type = "prob")$p[,1])
fitnm2<-mean(predict(ml2, newdata = mewdn, type = "prob")$p[,2])
fitym2 <- mean(predict(ml2, newdata = mewdy, type = "prob")$p[,2])
fitns2<-mean(predict(ml2, newdata = mewdn, type = "prob")$p[,3])
fitys2<- mean(predict(ml2, newdata = mewdy, type = "prob")$p[,3])
efm2<-exp(log(fitym2/fityc2)-log(fitnm2/fitnc2))
efs2<-exp(log(fitys2/fityc2)-log(fitns2/fitnc2))

mewdn <- ind[ind$nsfh==1,]
mewdn$unid<-0
mewdy <- ind[ind$nsfh==1,]
mewdy$unid<-1
fitnc3<-mean(predict(ml3, newdata = mewdn, type = "prob")$p[,1])
fityc3 <- mean(predict(ml3, newdata = mewdy, type = "prob")$p[,1])
fitnm3<-mean(predict(ml3, newdata = mewdn, type = "prob")$p[,2])
fitym3 <- mean(predict(ml3, newdata = mewdy, type = "prob")$p[,2])
fitns3<-mean(predict(ml3, newdata = mewdn, type = "prob")$p[,3])
fitys3<- mean(predict(ml3, newdata = mewdy, type = "prob")$p[,3])
efm3<-exp(log(fitym3/fityc3)-log(fitnm3/fitnc3))
efs3<-exp(log(fitys3/fityc3)-log(fitns3/fitnc3))

mewdn <- ind[ind$nsfh==0,]
mewdn$unid<-0
mewdy <- ind[ind$nsfh==0,]
mewdy$unid<-1
fitnc4<-mean(predict(ml4, newdata = mewdn, type = "prob")$p[,1])
fityc4 <- mean(predict(ml4, newdata = mewdy, type = "prob")$p[,1])
fitnm4<-mean(predict(ml4, newdata = mewdn, type = "prob")$p[,2])
fitym4 <- mean(predict(ml4, newdata = mewdy, type = "prob")$p[,2])
fitns4<-mean(predict(ml4, newdata = mewdn, type = "prob")$p[,3])
fitys4<- mean(predict(ml4, newdata = mewdy, type = "prob")$p[,3])
efm4<-exp(log(fitym4/fityc4)-log(fitnm4/fitnc4))
efs4<-exp(log(fitys4/fityc4)-log(fitns4/fitnc4))

#Put in a table
sink("mpc.Rnw")

cat(paste(
  '\\footnotesize 
\\begin{tabular}{@{\\extracolsep{5pt}}lcccc} 
\\\\[-1.8ex]\\hline 
\\hline 
\\\\[-1.8ex] & \\multicolumn{1}{c}{Full Sample} & \\multicolumn{1}{c}{Resident}& \\multicolumn{1}{c}{NSFH}& \\multicolumn{1}{c}{NSFG} \\\\ 
\\\\[-1.8ex] & \\multicolumn{1}{c}{(1)} & \\multicolumn{1}{c}{(2)} & \\multicolumn{1}{c}{(3)} & \\multicolumn{1}{c}{(4)}\\\\ 
\\hline \\\\[-1.8ex] 
\\\\[-2.2ex] & \\multicolumn{4}{c}{Risk of Marriage relative to Cohabitation} \\\\  
 \\hline \\\\[-1.8ex]
 Unilateral Divorce & ',round(p1a,digits=2),' & ',round(p2a,digits=2),' & ',round(p3a,digits=2),' & ',round(p4a,digits=2),' \\\\ 
  & (',round(s1a,digits=2),') & (',round(s2a,digits=2),') & (',round(s3a,digits=2),') & (',round(s4a,digits=2),') \\\\  
 \\hline \\\\[-1.8ex]
 Average Relative Risk & ',round(efm1,digits=2),' & ',round(efm2,digits=2),' & ',round(efm3,digits=2),' & ',round(efm4,digits=2),' \\\\ 
 \\hline \\\\[-1.8ex]
 \\\\[-2.2ex] & \\multicolumn{4}{c}{Risk of Separation relative to Cohabitation} \\\\  
 \\hline \\\\[-1.8ex]
 Unilateral Divorce & ',round(p1b,digits=2),' & ',round(p2b,digits=2),' & ',round(p3b,digits=2),' & ',round(p4b,digits=2),' \\\\ 
  & (',round(s1b,digits=2),') & (',round(s2b,digits=2),') & (',round(s3b,digits=2),') & (',round(s4b,digits=2),') \\\\  
 \\hline \\\\[-1.8ex]
 Average Relative Risk & ',round(efs1,digits=2),' & ',round(efs2,digits=2),' & ',round(efs3,digits=2),' & ',round(efs4,digits=2),' \\\\ 
 \\hline \\\\[-1.8ex]
State Fixed effects & Yes & Yes & Yes & Yes \\\\ 
Year Fixed effects & Yes & Yes & Yes & Yes \\\\ 
Age Polynomial & Yes & Yes & Yes & Yes \\\\
Picewise Duration & Yes & Yes & Yes & Yes \\\\ 
\\hline
Observations & \\multicolumn{1}{c}{',round(length(ind$id),digits=2),'} & \\multicolumn{1}{c}{',round(length(ind$id[ind$keep==1]),digits=2),'} & \\multicolumn{1}{c}{',round(length(ind$id[ind$nsfh==1]),digits=2),'} & \\multicolumn{1}{c}{',round(length(ind$id[ind$nsfh==0]),digits=2),'} \\\\ 
\\hline
Censored spells(\\%) & \\multicolumn{1}{c}{',round(100*mean(ind$censored[ind$sel==1]),digits=2),'} & \\multicolumn{1}{c}{',round(100*mean(ind$censored[ind$keep==1 & ind$sel==1]),digits=2),'} & \\multicolumn{1}{c}{',round(100*mean(ind$censored[ind$nsfh==1 & ind$sel==1]),digits=2),'} & \\multicolumn{1}{c}{',round(100*mean(ind$censored[ind$nsfh==0 & ind$sel==1]),digits=2),'} \\\\ 
\\hline 
\\hline \\\\[-1.8ex] 
\\end{tabular}'))

sink()
Sweave("mpc.Rnw")


#