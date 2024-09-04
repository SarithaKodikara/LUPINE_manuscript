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
load("Sim_Data/newSimData_n23.rds")
load("Sim_Data/newLibSize_n23.rds")
para<-readRDS("Simulation/Sim_para/para.rds")

n<-23
d<-54
iter<-50

clr_trans<-function(data){
  trans_data<-t(SpiecEasi::clr(t(data)))
  return(trans_data)
}

permutation_test<-function(residuals_i, residuals_j, iterations=100, seed=1234){
  perm_pvalues<-c()
  raw_cor<-cor(residuals_i, residuals_j)
  for(i in 1:iterations){
    perm_residuals_j <- sample(residuals_j)
    perm_pvalues[i] <- cor(residuals_i, perm_residuals_j)
  }
  perm_pvalue <- sum(raw_cor<perm_pvalues)/iterations
  
  return(perm_pvalue)
}



auc_res<-function( iter, day, method='sparcc'){
  #### Keeping track of the run
  cat(paste("\n","iter=", iter, "; ","day=",day,"; ","method=",method,"\n"), file="iteration.txt", append = TRUE)
  #### Getting the inferred network
  if(day<6){
    binary.graph<-para$binary.graph1_n
  }else{
    binary.graph<-para$binary.graph1_h
  }
  
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
        u1<-scale(X_1, center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<-lm(log(X_1[,taxa1[i]]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        r_j<-lm(log(X_1[,taxa2[i]]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        pcor[taxa1[i], taxa2[i]]<<-pcor[taxa2[i], taxa1[i]]<<-cor.test((r_i$residuals), 
                                                                     (r_j$residuals))$estimate
        pcor.pval[taxa1[i], taxa2[i]]<<-pcor.pval[taxa2[i], taxa1[i]]<<-1-cor.test((r_i$residuals), 
                                                                                 (r_j$residuals))$p.value
        
      })
    })
    res<-list(Estimate=pcor, pvalue=pcor.pval, RunTime=elap_time)
    auc_obj <- precrec::auc(evalmod(scores = pcor.pval[upper.tri(pcor.pval)], 
                                    labels = true_graph[upper.tri(true_graph)]))
  }else if(method=='PCA_log'){
    elap_time<-system.time({
      pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
      pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
      
      X_1<-log(newSimData[[iter]][,, day]+1)
      X_1_pca<-pca(X_1, ncomp=1, scale=TRUE)
      loadings_m<-matrix(c(X_1_pca$loadings$X),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(X_1, center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<-lm(X_1[,taxa1[i]]~u1,offset=log(newLibSize[[iter]][,day]))
        r_j<-lm(X_1[,taxa2[i]]~u1,offset=log(newLibSize[[iter]][,day]))
        pcor[taxa1[i], taxa2[i]]<<-pcor[taxa2[i], taxa1[i]]<<-cor.test((r_i$residuals), 
                                                                       (r_j$residuals))$estimate
        pcor.pval[taxa1[i], taxa2[i]]<<-pcor.pval[taxa2[i], taxa1[i]]<<-1-cor.test((r_i$residuals), 
                                                                                   (r_j$residuals))$p.value
        
      })
    })
    res<-list(Estimate=pcor, pvalue=pcor.pval, RunTime=elap_time)
    auc_obj <- precrec::auc(evalmod(scores = pcor.pval[upper.tri(pcor.pval)], 
                                    labels = true_graph[upper.tri(true_graph)]))
  }else if(method=='PCA_clr'){
    elap_time<-system.time({
      pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
      pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
      
      X_1<-clr_trans(newSimData[[iter]][,, day])
      X_1_pca<-pca(X_1, ncomp=1, scale=TRUE)
      loadings_m<-matrix(c(X_1_pca$loadings$X),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(X_1, center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<-lm(X_1[,taxa1[i]]~u1)
        r_j<-lm(X_1[,taxa2[i]]~u1)
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
      
      X_t <- newSimData[[iter]][,, (day)]
      
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
        res_netPLS1<-pls(d1$X1, X_t, 
                         ncomp=1, max.iter =1000)
        
      }else{
        res_netPLS1<-block.pls(d1, X_t, ncomp=1, max.iter =1000)
        
      }
      
      loadings_m<-matrix(c(res_netPLS1$loadings$Y),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(X_t, center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<<-lm(log(X_t[,taxa1[i]]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        r_j<<-lm(log(X_t[,taxa2[i]]+1)~u1,offset=log(newLibSize[[iter]][,day]))
        pcor[taxa1[i], taxa2[i]]<<-pcor[taxa2[i], taxa1[i]]<<-cor.test((r_i$residuals), 
                                                                     (r_j$residuals))$estimate
        pcor.pval[taxa1[i], taxa2[i]]<<-pcor.pval[taxa2[i], taxa1[i]]<<-1-cor.test((r_i$residuals), 
                                                                                 (r_j$residuals))$p.value
        
      })
    })
    
    res<-list( Estimate=pcor, pvalue=pcor.pval, RunTime=elap_time)
    auc_obj <- precrec::auc(evalmod(scores = pcor.pval[upper.tri(pcor.pval)], 
                                    labels = true_graph[upper.tri(true_graph)]))
  }else if(method=='bPLS_log'){
    
    elap_time<-system.time({
      pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
      pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
      
      X_t<-log(newSimData[[iter]][,, (day)]+1)
      
      day_i=day
      d1<-list(X1=log(newSimData[[iter]][,, (day_i-1)]+1))
      day_i=day_i-1
      index=3
      while(day_i!=1){
        d1<-list.append(d1,log(newSimData[[iter]][,, (day_i-1)]+1))
        names(d1)[index-1]<-paste("X", index-1, sep = "")
        index=index+1
        day_i=day_i-1 
      }
      if(day==2){
        res_netPLS1<-pls(d1$X1, X_t, 
                         ncomp=1, max.iter =1000)
        
      }else{
        res_netPLS1<-block.pls(d1, X_t, 
                               ncomp=1, max.iter =1000)
        
      }
      
      loadings_m<-matrix(c(res_netPLS1$loadings$Y),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(X_t, center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<<-lm(X_t[,taxa1[i]]~u1,offset=log(newLibSize[[iter]][,day]))
        r_j<<-lm(X_t[,taxa2[i]]~u1,offset=log(newLibSize[[iter]][,day]))
        pcor[taxa1[i], taxa2[i]]<<-pcor[taxa2[i], taxa1[i]]<<-cor.test((r_i$residuals), 
                                                                       (r_j$residuals))$estimate
        pcor.pval[taxa1[i], taxa2[i]]<<-pcor.pval[taxa2[i], taxa1[i]]<<-1-cor.test((r_i$residuals), 
                                                                                   (r_j$residuals))$p.value
        
      })
    })
    
    res<-list( Estimate=pcor, pvalue=pcor.pval, RunTime=elap_time)
    auc_obj <- precrec::auc(evalmod(scores = pcor.pval[upper.tri(pcor.pval)], 
                                    labels = true_graph[upper.tri(true_graph)]))
  }else if(method=='bPLS_clr'){
    
    elap_time<-system.time({
      pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
      pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
      
      X_t<-clr_trans(newSimData[[iter]][,, (day)])
      
      day_i=day
      d1<-list(X1=clr_trans(newSimData[[iter]][,, (day_i-1)]))
      day_i=day_i-1
      index=3
      while(day_i!=1){
        d1<-list.append(d1,clr_trans(newSimData[[iter]][,, (day_i-1)]))
        names(d1)[index-1]<-paste("X", index-1, sep = "")
        index=index+1
        day_i=day_i-1 
      }
      if(day==2){
        res_netPLS1<-pls(d1$X1, X_t, 
                         ncomp=1, max.iter =1000)
        
      }else{
        res_netPLS1<-block.pls(d1, X_t, 
                               ncomp=1, max.iter =1000)
        
      }
      
      loadings_m<-matrix(c(res_netPLS1$loadings$Y),ncol = 1)
      
      sapply(1:len, function(i){
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(X_t, center=TRUE, scale=TRUE)%*%loading_tmp
        r_i<<-lm(X_t[,taxa1[i]]~u1)
        r_j<<-lm(X_t[,taxa2[i]]~u1)
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
