
setwd("C:/Users/macristina.perez/Documents/UMassD/Classes/2022/fall_2022/length-based/Class 8/SC3")

library(dplyr)
library(ggplot2)
graphics.off()
rm(list=ls(all=TRUE))

#shell("SBFT.exe")

# model without selectivity and by estimating proportions at age, 
# mean length of the first age class, mean length of the last age class, 
# lambda 1 (magnitude of sds for mean length-at-age), 
# and lambda 2 (length-dependent trend in sd for mean length-at-age).


# How similar are your results to Fournier’s in Table 6? Do these results seem reasonable? 
# IOW (in other words), if you didn’t have another set of estimates to compare to, would you think these were 
# reasonable? Consult Table 7 for previously published values of K and Lin for this species.

# K is bigger (growth rate) and Linf is smaller 
# mean length age in the first age and mean lenght in the last age 
#Visualize results
MFLFs=read.delim("sbft.rep",sep="",header=FALSE)
otherres=read.delim("other.rep",sep="",header=FALSE)
obsLF=otherres[1:10,]
obsNAL=otherres[11:20,]
#Prep for visualation
midbins=seq(30.5,208.5,2)
ts=10 #number of datasets
obsLFplot=t(rbind(obsLF,midbins,rep(2,length(midbins))))
colnames(obsLFplot)=c(seq(1,ts),"Midbin","Model")
MFLFsplot=t(rbind(MFLFs,midbins,rep(3,length(midbins))))
colnames(MFLFsplot)=c(seq(1,ts),"Midbin","Model")
resLF=data.frame(rbind(MFLFsplot,obsLFplot))
resLF$Model=as.factor(resLF$Model);levels(resLF$Model)= c("Obs","Est")
#Plot observed vs estimated LFs
graphics.off()
for(t in 1:ts){
  resLFsub=resLF[,t]
  resLFsub=data.frame(cbind(resLFsub,resLF[,ts+1],resLF[,ts+2]))
  colnames(resLFsub)=c("Proportion","Midbin","Model")
  resLFsub$Model=as.factor(resLFsub$Model)
  levels(resLFsub$Model) = c("Obs","Est")
  g=ggplot(resLFsub,aes(x=Midbin,y=Proportion,color=Model))+
    geom_line()+
    geom_point(size=4,aes(shape=Model))+
    scale_shape_manual(values=c(16,15))+
    scale_color_manual(values=c("#FF33B5","#660882"))+
    theme_bw(base_size = 20)+
    xlab("Length(midbin)")+
    theme(
      legend.position = c(.95, .95),
      legend.justification = c("right", "top"),
      legend.box.just = "right",
      legend.margin = margin(6, 6, 6, 6)
    )
  print(g)
}
#Grab vonB parameter estimates
vbkest=otherres[21,1]
linf=otherres[22,1]
mLfirst=otherres[23,1]
mLlast=otherres[24,1]
vbkest
linf
mLfirst
mLlast

# How sensitive is the model to the assumption of number of ages? How sensitive is the model to starting
# values?

# con 10 edades lINF  172   Y k 17
# CON 20 EDADES LINF 481 Y K 0.014
# con 5 edades Linf 2797 y k 0.006

# con prior k de 0.17  LINF 172.749 Y K 0.17604
# con prior k de 0.25  LINF 172.734 Y K 0.1757

# con prior k de 0.25 Y PRIOR ALFA 3  LINF 172.736 Y K 0.175631