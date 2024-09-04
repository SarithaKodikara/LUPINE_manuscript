library(MASS)
library(SpiecEasi)
library(corpcor)
library(tseries)
library(dplyr)
library(ggplot2)
library(mixOmics)
library(tidyr)
library(faux)
library(gplots)
library(simstudy)
library(fake)
library(MGLM)
library(lcmix)


multi.INAR1.sim<-function(d, t, meanV1, arPara1, Cor1){
  # Initialize variables
  X <- matrix(NA, nrow = t, ncol = d)
  
  rateInnov1<-meanV1
  
  X[1,1:d] <- genCorGen_Ex(1, nvars =d, params1 = rateInnov1, dist = "poisson", 
                           corMatrix =Cor1 )$X  # Initial value
  
  for (i in 2:(t)) {
      X[i, 1:d] <- genCorGen_Ex(1, nvars =d, params1 = rateInnov1, 
                                dist = "poisson", corMatrix =Cor1 )$X + 
        genCorGen_Ex(1, nvars =d, params1 = X[(i-1),1:d], 
                     params2 = rep(arPara1,d), dist = "binomial", 
                     corMatrix =Cor1)$X
    
  }
  return(X)
  
}

# n:= number of Individuals
# d:= number of OTUs
# t:= number of Timepoints
# rho:= parameter used in AR(1)
# seed:= random seed to generate the network
SimCount<-function(n, d, t, seed, cor1, gamma_par1, arPara1=0.5){
  
  set.seed(seed)
  
  # setup vector for data
  Y.sim <- array(NA, c(n, d, t))
  
  # Initial mean OTU values generated from a multivariate gamma for each individuals
  shape1<-para$gamma_par_n[,1]
  scale1<-para$gamma_par_n[,2]
  rate1<-1/scale1
  mu_subject1<-rmvgamma(n, shape=shape1, rate=rate1, corr=cor1)
  
  for(i in 1:n){
    Y.sim[i,1:d,1:t]<-t(multi.INAR1.sim(d=d, t=t, 
                                        meanV1=mu_subject1[i,], arPara1, cor1))
  }
  
  
  
  return(Y.sim)
}

