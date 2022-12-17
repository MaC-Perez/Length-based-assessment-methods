setwd("C:/Users/macristina.perez/Documents/UMassD/Classes/2022/fall_2022/length-based/Class 2")

library(TropFishR)
library(fishmethods)
library(MLZ)
library(TMB)
library(dplyr)
graphics.off()
rm(list=ls(all=TRUE))


###################
### EXCERSISE 1 ###
###################

#calculate vB growth and BH length based catch curve (Yellow Pearch example)

YPdat=read.csv("YPlens.csv")
colnames(YPdat)=c("age","length")
growmod=growth(intype=1,unit=1,size=YPdat$length,age=YPdat$age,
               calctype=1,wgtby=1,error=1,Sinf=300,K=0.3,t0=-1)
growmod$vout

#Nonlinear regression model
#model: size ~ Sinf * (1 - exp(-(K * (age - t0))))
#data: x
#Sinf        K       t0 
#315.8946   0.2678  -1.5118 
#residual sum-of-squares: 2503805

#Number of iterations to convergence: 3 
#Achieved convergence tolerance: 3.543e-08


#compare with the results from the excel file 

###################
### EXCERSISE 2 ###
###################

data(synLFQ2)
synLFQ2

Z1960=Z_BevertonHolt(synLFQ2, catch_columns = 1, Lprime_tprime = 47.5)
Z1970=Z_BevertonHolt(synLFQ2, catch_columns = 2, Lprime_tprime = 47.5)
Z1980=Z_BevertonHolt(synLFQ2, catch_columns = 3, Lprime_tprime = 47.5)
Z1960$Z
Z1970$Z
Z1980$Z
BHests=c(Z1960$Z,Z1970$Z,Z1980$Z)
BHests



###################
### EXCERSISE 3 ###
###################


synLFQ2mod=synLFQ2[-c(1,5)]
synLFQ2mod$t0=0
synLFQ2mod$catch=synLFQ2$catch[]
dev.new()
PW1960=powell_wetherall(synLFQ2mod,catch_columns=1)
PW1970=powell_wetherall(synLFQ2mod,catch_columns=2)
PW1980=powell_wetherall(synLFQ2mod,catch_columns=3)
PWests=c(PW1960$ZK*PW1960$K,PW1970$ZK*PW1970$K,PW1980$ZK*PW1980$K)


###################
### EXCERSISE 4 ###
###################


LC1960 <- catchCurve(synLFQ2mod,catch_columns=1)
LC1970 <- catchCurve(synLFQ2mod,catch_columns=2)
LC1980 <- catchCurve(synLFQ2mod,catch_columns=3)
LCests=c(LC1960$Z,LC1970$Z,LC1980$Z)


###################
### EXCERSISE 5 ###
###################


yrs=c(1960,1970,1980)
allres=cbind(yrs,BHests,PWests,LCests)
plot(yrs,PWests,type="b",col="blue",ylim=c(0.5,1.5),ylab="Z or relative Z")
lines(yrs,BHests,type="b",col="red",pch="*")
lines(yrs,LCests,type="b",col="black",pch="+")
legend(1960,1.5,legend=c("PW","BH","LC"),
       col=c("blue","red","black"),pch=c("o","*","+"), ncol=1)


###################
### EXCERSISE 6 ###
###################

data(SilkSnapper)
SS.dataset= new("MLZ_data", Year = 1983:2013,
                Len_df = SilkSnapper, length.units = "mm",vbLinf = 794, vbK = 0.1)
plot(SS.dataset, type = "comp")
