#get the final sample for the analysis


## remove (almost) everything in the working environment.
rm(list = ls())


list.of.packages <- c("haven", "mcmc","coda","MNP","stargazer","foreing","nnet","ggplot2","reshape2")

ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE, repos='http://cran.us.r-project.org')
sapply(pkg, require, character.only = TRUE)
}

ipak(list.of.packages)


#for reading .dta file

library(haven)
library(readstata13)
library(mcmc)
library(coda)

#for multinomial logit+tobit+probit
#library(mlogit)
#library(mnlogit)
#library(margins)
#library(gmnl)
#library(VGAM)
library(MNP)
#library(sampleSelection)
#library(survival)
#library(rms)

#fixed effects regression
#library(plm)
#require(lmtest)

#for exporting in Latex
library(stargazer)
#require(foreign)
#require(nnet)
#require(ggplot2)
#require(reshape2)


##################################################################
##################################################################
####  SAMPLE SELECTION AS A FIRST THING
##################################################################
##################################################################
directory<-"C:\\Users\\Fabio\\Dropbox\\JMP\\presentation\\phd_apero_08_2019"
#directory<-"C:\\Users\\Fabio\\Dropbox\\cohabitation and divorce laws\\coh_text"
setwd(directory)


##################################################################
##################################################################
####  COHABITATION LENGTH
##################################################################
##################################################################

#ind <-read_dta("C:/Users/Fabio/Dropbox/JMP/empirical analysis/NSFG/88/NSFG88s2.dta")

ind <-read.dta13("C:/Users/Fabio/Dropbox/JMP/empirical analysis/NSFG/88/NSFG88s2.dta")




ml1 <- mnp(rel ~  unid+factor(st)+factor(timey),
           data =ind, n.draws = 50, base="2", verbose = TRUE)

ml2 <- mnp(rel ~  unid+factor(st)+factor(timey),
           data =ind, n.draws = 50,base="2",
           coef.start=rep(1,254),verbose = TRUE)

ml3 <- mnp(rel ~  unid+factor(st)+factor(timey),
           data =ind, n.draws = 50,base="2",  verbose = TRUE)


ml4 <- mnp(rel ~  unid+factor(st)+factor(timey),
           data =ind, n.draws = 50,base="2",  verbose = TRUE)


res.coda <- mcmc.list(chain1=mcmc(ml1$param),
                      chain2=mcmc(ml2$param))
                                  



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
beastn<-predict(ml1, newdata = mewdn, type = "prob")
beasty<-predict(ml1, newdata = mewdn, type = "prob")

fitnc1<-mean(beastn$p[,1])
fityc1 <- mean(beasty$p[,1])
fitnm1<-mean(beastn$p[,2])
fitym1 <- mean(beasty$p[,2])
fitns1<-mean(beastn$p[,3])
fitys1<- mean(beasty$p[,3])
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
sink("mpcc.Rnw")

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
\\hline \\\\[-1.8ex] 
\\end{tabular}'))

sink()
Sweave("mpcc.Rnw")

