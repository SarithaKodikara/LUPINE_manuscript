library(e1071)
library(ade4)
library(graphsim)
library(vegan)

load(paste0("Results/PPT/LUPINE_Day2.rdata"))
net1<-(res$pvalue<0.05)*1
net1<-apply(net1,c(1,2), function(x){ifelse(is.na(x),0,x)})

load(paste0("Results/MDE/LUPINE_Day2.rdata"))
net2<-(res$pvalue<0.05)*1
net2<-apply(net2,c(1,2), function(x){ifelse(is.na(x),0,x)})


lapl1<-as.dist(hamming.distance(net1))
lapl2<-as.dist(hamming.distance(net2))
mantel.rtest(lapl1,lapl2)

