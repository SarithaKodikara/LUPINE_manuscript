library(abind)
library(mixOmics)
library(SpiecEasi)
library(Matrix)
library(igraph)
library(corpcor)
library(mixOmics)
library(caret)
library(pls)
library(lars)
library(glassoFast)
library(patchwork)
library(ggplotify)
library(dplyr)
library(Hmisc)
library(rlist)
library(precrec)

load(paste0("Sim_Data/newSimData_n23_e",Edge_percentage,".rds"))
load(paste0("Sim_Data/newLibSize_n23_e",Edge_percentage,".rds"))
load(paste0("Sim_Data/graph",Edge_percentage,".rds"))

n<-23
d<-54
iter<-50


auc_res<-function( iter, day, method='sparcc'){
  #### Keeping track of the run
  cat(paste("\n","iter=", iter, "; ","day=",day,"; ","method=",method,"\n"), file="iteration.txt", append = TRUE)
  #### Getting the inferred network
  binary.graph<-graph_new
  
  true_graph<-structure(binary.graph,class="graph")
  
  ig.org <- adj2igraph(true_graph)%>%
    set_vertex_attr("name", value = paste0(1: dim(binary.graph)[1]))
  
  Y_count<-newSimData[[iter]][,,day]
  nOTU<-dim(Y_count)[2]
  
  len<-nOTU*(nOTU-1)/2
  taxa1=unlist(lapply(1:nOTU, function(i)rep(i,(nOTU-i))))
  taxa2=unlist(lapply(2:nOTU, function(i)seq(i,nOTU,1)))
  
  if(method=='PCA'){
    elap_time<-system.time({
      pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
      pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
      
      X_1<-newSimData[[iter]][,, day]
      X_1_pca<-pca(X_1, ncomp=1, scale=TRUE)
      loadings_m<-matrix(c(X_1_pca$loadings$X),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(newSimData[[iter]][,, day], center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<-lm(log(newSimData[[iter]][,taxa1[i], day]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        r_j<-lm(log(newSimData[[iter]][,taxa2[i], day]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        pcor[taxa1[i], taxa2[i]]<<-pcor[taxa2[i], taxa1[i]]<<-cor.test((r_i$residuals), 
                                                                     (r_j$residuals))$estimate
        pcor.pval[taxa1[i], taxa2[i]]<<-pcor.pval[taxa2[i], taxa1[i]]<<-1-cor.test((r_i$residuals), 
                                                                                 (r_j$residuals))$p.value
        
      })
    })
    res<-list(Estimate=pcor, pvalue=pcor.pval, RunTime=elap_time)
    auc_obj <- precrec::auc(evalmod(scores = pcor.pval[upper.tri(pcor.pval)], 
                                    labels = true_graph[upper.tri(true_graph)]))
  }else if(method=='bPLS'){
    
    elap_time<-system.time({
      pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
      pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
      
      day_i=day
      d1<-list(X1=newSimData[[iter]][,, (day_i-1)])
      day_i=day_i-1
      index=3
      while(day_i!=1){
        d1<-list.append(d1,newSimData[[iter]][,, (day_i-1)])
        names(d1)[index-1]<-paste("X", index-1, sep = "")
        index=index+1
        day_i=day_i-1 
      }
      if(day==2){
        res_netPLS1<-pls(d1$X1, newSimData[[iter]][,,day], 
                         ncomp=1, max.iter =1000)
        
      }else{
        res_netPLS1<-block.pls(d1, newSimData[[iter]][,,day], ncomp=1, max.iter =1000)
        
      }
      
      loadings_m<-matrix(c(res_netPLS1$loadings$Y),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(newSimData[[iter]][,, day], center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<<-lm(log(newSimData[[iter]][,taxa1[i], day]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        r_j<<-lm(log(newSimData[[iter]][,taxa2[i], day]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        pcor[taxa1[i], taxa2[i]]<<-pcor[taxa2[i], taxa1[i]]<<-cor.test((r_i$residuals), 
                                                                     (r_j$residuals))$estimate
        pcor.pval[taxa1[i], taxa2[i]]<<-pcor.pval[taxa2[i], taxa1[i]]<<-1-cor.test((r_i$residuals), 
                                                                                 (r_j$residuals))$p.value
        
      })
    })
    
    res<-list( Estimate=pcor, pvalue=pcor.pval, RunTime=elap_time)
    auc_obj <- precrec::auc(evalmod(scores = pcor.pval[upper.tri(pcor.pval)], 
                                    labels = true_graph[upper.tri(true_graph)]))
    }
  
  roc_prc<-data.frame(Iter=iter, Day=day, Method=method, ROC=auc_obj[1,4], PRC=auc_obj[2,4])
  return(roc_prc)
  
}
