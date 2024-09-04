library(RColorBrewer)
library(patchwork)
library(igraph)
library(ggplotify)
library(abind)
library(qgraph)
library(circlize)
library(ComplexHeatmap)
library(vegan)
#### Reading Data ####
source("utilities.R")
source('SimFunction.R') 
para<-readRDS("Sim_para/para.rds")
Colors=rev(brewer.pal(11,"Spectral"))
Colors=colorRampPalette(Colors)(100)

Edge_percentage = 0.05


getCorPcor<-function(graph_i, seed){
  set.seed(seed)
  graph<-structure(graph_i, class='graph')
  pre<-graph2prec(graph)
  library(qgraph)
  pcor<-as.matrix(wi2net(pre))
  (min(pcor));(max(pcor))
  diag(pcor)<-1
  test1<-is.positive.definite(pcor)
  test2<-isSymmetric(pcor)
  
  cor<-cov2cor(prec2cov(pre))
  test3<-is.positive.definite(cor)
  test4<-isSymmetric(cor)
  
  while(!(test1&test2&test3&test4)){
    graph<-structure(graph_i, class='graph')
    pre<-graph2prec(graph)
    library(qgraph)
    pcor<-as.matrix(wi2net(pre))
    (min(pcor));(max(pcor))
    diag(pcor)<-1
    test1<-is.positive.definite(pcor)
    test2<-isSymmetric(pcor)
    
    cor<-cov2cor(prec2cov(pre))
    test3<-is.positive.definite(cor)
    test4<-isSymmetric(cor)
  }
  return(list(cor=cor, pcor=pcor))
}

modifyGraph<-function(graph, percen=0.005){
  numberofEdges<- sum(graph)/2
  newEdges <- round((dim(graph)[1]*(dim(graph)[1]-1)/2)*percen)
  newgraph <- graph
  if(numberofEdges>newEdges){
    con_index<-which(newgraph*upper.tri(newgraph)==1)
    index_to_0<-sample(con_index, numberofEdges-newEdges)
    newgraph[index_to_0] <- 0 
    newgraph[lower.tri(newgraph)] <-t(newgraph)[lower.tri(newgraph)]
  }else if(numberofEdges<newEdges){
    con_index<-which(newgraph*upper.tri(newgraph)==0)
    index_to_1<-sample(con_index, newEdges-numberofEdges)
    newgraph[index_to_1] <- 1 
    newgraph[lower.tri(newgraph)] <-t(newgraph)[lower.tri(newgraph)]
  }
  return(newgraph)
}

#Function that simulates counts nSim times and that have no NANs
set.seed(156431)
nSim=50
n=23
d=54
t=5
gamma_par1=para$gamma_par_n
graph1=para$binary.graph1_n
arPara1=0.5
 
graph_new<-modifyGraph(graph1, Edge_percentage)
save(graph_new,file = paste0("Sim_Data/graph",Edge_percentage,".rds"))

ig.org <- adj2igraph(graph_new)%>%
  set_vertex_attr("name", value = paste0(1:d))

png(paste0("Sim_Data/Plots/net_",Edge_percentage,".png"), width=1200, height=1200, res=120) # start export
plot.igraph(ig.org,vertex.size=6, edge.width=3,vertex.color="orange")
dev.off()

seed=sample(1:1000000,nSim)

nSim_count<-function(s){
  print(s)
  normal_cor_pcor<-getCorPcor(graph_new, seed[s])
  cor1<-normal_cor_pcor$cor
  pcor1<-normal_cor_pcor$pcor
  simData_i<-SimCount(n, d, t, seed[s], cor1, gamma_par1, arPara1)
  return(list(cor=cor1, pcor=pcor1, simData_i=simData_i))
}

allData<-lapply(1:nSim,function(s){nSim_count(s)})
SimData<-lapply(1:nSim,function(s){allData[[s]]$simData_i})
corData<-lapply(1:nSim,function(s){allData[[s]]$cor})
pcorData<-lapply(1:nSim,function(s){allData[[s]]$pcor})

save(SimData,file = paste0("Sim_Data/SimData_n23_e",Edge_percentage,".rds"))
save(corData,file = paste0("Sim_Data/cor_e",Edge_percentage,".rds"))
save(pcorData,file = paste0("Sim_Data/pcor_e",Edge_percentage,".rds"))

LibSize<-lapply(1:nSim,function(s){apply(SimData[[s]],c(1,3),sum)})

save(LibSize,file = paste0("Sim_Data/LibSize_n23_e",Edge_percentage,".rds"))

set.seed(1234)
newLibSize<-lapply(1:nSim,function(s){apply(LibSize[[s]],c(1,2), function(x){sample(5000:8000, 1)})})

save(newLibSize,file = paste0("Sim_Data/newLibSize_n23_e",Edge_percentage,".rds"))
newSimData<-lapply(1:nSim,function(s){array(unlist(lapply(1:t, 
                            function(day){rrarefy(x=SimData[[s]][,,day],sample=newLibSize[[s]][,day])})),dim=c(n,d,t))})
save(newSimData,file = paste0("Sim_Data/newSimData_n23_e",Edge_percentage,".rds"))
