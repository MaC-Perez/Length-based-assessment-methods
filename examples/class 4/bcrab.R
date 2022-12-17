setwd("C:/Users/macristina.perez/Documents/UMassD/Classes/2022/fall_2022/length-based/Class 4")
source("read.report.R")

library(ggplot2)
library(reshape2)
graphics.off()
rm(list=ls(all=TRUE))

shell("bcrab.exe")

file.rename("bcrab.rep","bcrab.txt")
repfile=read.delim("bcrab.txt",header=FALSE,sep=" ")
popmod_ests=repfile[1:31,]
colnames(popmod_ests)=c("Year","Adult_N","Rec_N","u","log_R_devs")
obsmod_ests=repfile[32:62,]
colnames(obsmod_ests)=c("Year","Est_adult","Obs_adult","Est_juv","Obs_juv")
popmod_estsN=melt(popmod_ests[,1:3],id=c("Year"))
compad=melt(obsmod_ests[,1:3],id="Year")
colnames(compad)=c("Year","ObsPred","Est")
compad$ObsPred=factor(compad$ObsPred)
compjuv=melt(obsmod_ests[,c(1,4:5)],id="Year")
colnames(compjuv)=c("Year","ObsPred","Est")
compjuv$ObsPred=factor(compjuv$ObsPred)
catch=read.delim("bcrab_catch.txt",header=T)

###################################
ggplot(data=popmod_ests, aes(x=Year, y=Adult_N)) +
  geom_line()+
  geom_point()+
  labs(x="Year", y = "Adult abundance")+
  theme_classic()

ggplot(data=popmod_ests, aes(x=Year, y=Rec_N)) +
  geom_line()+
  geom_point()+
  labs(x="Year", y = "Recruit abundance")+
  theme_classic()

ggplot(data=popmod_ests, aes(x=Year, y=log_R_devs)) +
  geom_line()+
  geom_point()+
  labs(x="Year", y = "Recruitment deviations",ylim=c(-0.6,0.6))+
  geom_hline(yintercept=0,col="red")+
  theme_classic()

##################################################
ggplot(data=popmod_ests, aes(x=Year, y=u)) +
  geom_line()+
  geom_point()+labs(x="Year", y = "Exploitation Rate")+
  theme_classic()

ggplot(data=catch, aes(x=Year, y=Catch)) +
  geom_line()+
  geom_point()+
  labs(x="Year", y = "Catch")+
  theme_classic()
#####################################################

ggplot(data=compad, aes(x=Year, y=Est,group=ObsPred)) +
  geom_line(aes(col=ObsPred))+
  geom_point(aes(shape=ObsPred))+
  labs(x="Year", y = "Adult index")+
  theme_classic()

ggplot(data=compjuv, aes(x=Year, y=Est,group=ObsPred)) +
  geom_line(aes(col=ObsPred))+
  geom_point(aes(shape=ObsPred))+
  labs(x="Year", y = "Recruit index")+
  theme_classic()

###################################################################

std<-read.table("bcrab.std",skip=-1)
std<-std[-1,]
names(std)<-c("index","name","value","std.dev")

log_R_devs  <-which(std$name=="log_R_devs")
log_R_devs  <-std[log_R_devs,]
log_R_devs  <-as.numeric(log_R_devs  [,3])

rr<-dim(log_R_devs)
rpr.c<-as.numeric(levels(log_R_devs[,3]))[as.numeric(log_R_devs[,3])]
sup<-rpr.c+2*as.numeric(levels(log_R_devs[,4]))[as.integer(log_R_devs[,4])]
inf<-rpr.c-2*as.numeric(levels(log_R_devs[,4]))[as.integer(log_R_devs[,4])]
poli<-rbind(cbind(year,sup),cbind(rev(year),rev(inf)))


