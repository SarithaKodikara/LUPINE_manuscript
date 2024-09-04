library(RColorBrewer)
library(patchwork)
library(igraph)
library(ggplotify)
library(abind)
library(qgraph)
library(circlize)
library(ComplexHeatmap)
library(vegan)
library(parallel)
nworkers <- parallel::detectCores()
#### Reading Data ####
set.seed(156431)
source("utilities.R")
source('SimFunction.R') 
load("para_networks/graph.rds")
load("para_networks/cor_pcor_list.rds")
para<-readRDS("Sim_para/para.rds")
nSim=50
d=54
t=5
gamma_par1=para$gamma_par_n
arPara1=0.5
n=5000
seed=sample(1:1000000,nSim)

SimData<-mclapply(1:nSim, function(s){
  cat(paste("\n","size=", n, "; ","iter=",s,"\n"), file=paste0("iteration_sim",n,".txt"), append = TRUE)
  SimCount(n, d, t, seed[s], cor_pcor_list[[s]]$cor, gamma_par1, arPara1)}, 
  mc.cores = nworkers-2)

save(SimData,file = paste0("Sim_Data/SimData_n",n,".rds"))

LibSize<-lapply(1:nSim,function(s){apply(SimData[[s]],c(1,3),sum)})
save(LibSize,file = paste0("Sim_Data/LibSize_n",n,".rds"))

set.seed(1234)
newLibSize<-lapply(1:nSim,function(s){apply(LibSize[[s]],c(1,2), function(x){sample(5000:8000, 1)})})

save(newLibSize,file = paste0("Sim_Data/newLibSize_n",n,".rds"))

newSimData<-lapply(1:nSim,function(s){array(unlist(lapply(1:t, 
                                                          function(day){rrarefy(x=SimData[[s]][,,day],sample=newLibSize[[s]][,day])})),dim=c(n,d,t))})
save(newSimData,file = paste0("Sim_Data/newSimData_n",n,".rds"))

