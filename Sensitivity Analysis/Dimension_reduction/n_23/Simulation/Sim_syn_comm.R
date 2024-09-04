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


#####Visualisation of Precision & Correlation
Sim_plots<-function(Graph_org, Cor, pCor, iteration, stage){
  nOTU<- ncol(Cor)
  
  ig.org <- adj2igraph(Graph_org)%>%
    set_vertex_attr("name", value = paste0(1:nOTU))
  iso.org <- V(ig.org)[igraph::degree(ig.org)==0]
  g.org <- delete_vertices(ig.org, iso.org)
  
  png(paste("Sim_Data/Plots/Sim_", iteration,"_s",stage,"_p1.png", sep = ""), width=1200, height=1200, res=120) # start export
  plot.igraph(ig.org,vertex.size=6, edge.width=3,vertex.color="orange")
  dev.off()
  
  con<-data.frame(partial_Correlation=unlist(lapply(
    lapply(strsplit(as_ids(E(g.org)), "[|]"), as.numeric), 
    function(i)pCor[i[1],i[2]])))
  p1<-ggplot(con, aes(x= partial_Correlation)) +
    geom_histogram(data = con, fill = "green", alpha = 0.5, binwidth = 0.1)+labs(title="Connections and non-connections")
  p11<-as.ggplot(ComplexHeatmap::Heatmap(pCor, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(-1, 0, 1), c("#009E73", "white", "#D55E00")) ))
 
  
  
  ig.inf <- adj2igraph(abs(Cor)>0)%>%
    set_vertex_attr("name", value = paste0(1:nOTU))
  iso.inf <- V(ig.inf)[igraph::degree(ig.inf)==0]
  g.inf <- delete.vertices(ig.inf, iso.inf)
  
  TrueCor_FP<-data.frame(Correlation=unlist(lapply(
    lapply(strsplit(as_ids(E(g.inf)[!(as_ids(E(g.inf))%in%as_ids(E(g.org)))]), "[|]"), as.numeric), 
    function(i)Cor[i[1],i[2]])))
  TrueCor_TP<-data.frame(Correlation=unlist(lapply(
    lapply(strsplit(as_ids(E(g.inf)[as_ids(E(g.inf))%in%as_ids(E(g.org))]), "[|]"), as.numeric), 
    function(i)Cor[i[1],i[2]])))
  CombinedDataTrue <- rbind(TrueCor_FP, TrueCor_TP)
  p2<-ggplot(CombinedDataTrue, aes(x= Correlation)) + geom_histogram(data = TrueCor_FP, fill = "purple", alpha = 0.2, binwidth = 0.1) + 
    geom_histogram(data = TrueCor_TP, fill = "green", alpha = 0.5, binwidth = 0.1)+labs(title="Connections and non-connections")
   p21<-as.ggplot(ComplexHeatmap::Heatmap(Cor, cluster_rows = F, cluster_columns = F, col = colorRamp2(c(-1, 0, 1), c("#009E73", "white", "#D55E00")) ))
  
  
  
  png(paste("Sim_Data/Plots/Sim_",iteration,"_s",stage,"_p2.png", sep = ""), width=1200, height=1200, res=120) # start export
  print((p1|p11)/(p2|p21))
  dev.off()
}

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

#Function that simulates counts nSim times and that have no NANs
set.seed(156431)
nSim=50
n=23
d=54
t=10
interven_t= 6
gamma_par1=para$gamma_par_n
gamma_par2=para$gamma_par_h
graph1=para$binary.graph1_n
graph2=para$binary.graph1_h
arPara1=0.5
arPara2=0.1


seed=sample(1:1000000,nSim)

nSim_count<-function(s){
  print(s)
  normal_cor_pcor<-getCorPcor(graph1, seed[s])
  cor1<-normal_cor_pcor$cor
  pcor1<-normal_cor_pcor$pcor
  Sim_plots(Graph_org=graph1, Cor=cor1, pCor=pcor1, iteration=s, stage=1)
  
  hfhs_cor_pcor<-getCorPcor(graph2, seed[s])
  cor2<-hfhs_cor_pcor$cor
  pcor2<-hfhs_cor_pcor$pcor
  Sim_plots(Graph_org=graph2, Cor=cor2, pCor=pcor2, iteration=s, stage=2)
  
  simData_i<-SimCount(n, d, t, seed[s], cor1, gamma_par1, arPara1,
                      cor2, gamma_par2, arPara2, interven_t)
  return(simData_i)
}

SimData<-lapply(1:nSim,function(s){nSim_count(s)})


save(SimData,file = "Sim_Data/SimData_n23.rds")

LibSize<-lapply(1:nSim,function(s){apply(SimData[[s]],c(1,3),sum)})

save(LibSize,file = "Sim_Data/LibSize_n23.rds")

set.seed(1234)
newLibSize<-lapply(1:nSim,function(s){apply(LibSize[[s]],c(1,2), function(x){sample(5000:8000, 1)})})

save(newLibSize,file = "Sim_Data/newLibSize_n23.rds")
newSimData<-lapply(1:nSim,function(s){array(unlist(lapply(1:t, 
                            function(day){rrarefy(x=SimData[[s]][,,day],sample=newLibSize[[s]][,day])})),dim=c(n,d,t))})
save(newSimData,file = "Sim_Data/newSimData_n23.rds")
