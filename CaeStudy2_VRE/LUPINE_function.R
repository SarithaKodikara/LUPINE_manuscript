library(mixOmics)
library(dplyr)
library(rlist)

clr_trans<-function(data){
  trans_data<-t(SpiecEasi::clr(t(data)))
  return(trans_data)
}


## n:-number of participants; p:- number of taxa; t:- number of time points)
# data: an array of size n x p x t
# day: an numeric value [time point when the network is inferred]
# excluded_taxa: list of size t [taxa to be excluded when inferring the network at each time point]
# transform: TRUE [do a clr transformation on the data]
# lib_size: matrix of size n x t
# single: FALSE [longitudinal version of LUPINE]

LUPINE<- function(data, day, excluded_taxa=NULL, transform=FALSE, lib_size=NULL, single=FALSE){
  # Total number of taxa (p)
  nOTU_total <- dim(data)[2]
  # Extract matrix for current time point
  count_day <- data[,,day]
  # Extract taxa names after excluded taxa
  taxa_names<-colnames(count_day)[!(colnames(count_day)%in%excluded_taxa[[day]])]
  # pairwise taxa combinations
  nOTU<-length(taxa_names)
  len<-nOTU*(nOTU-1)/2
  taxa1=unlist(lapply(1:nOTU, function(i)rep(i,(nOTU-i))))
  taxa2=unlist(lapply(2:nOTU, function(i)seq(i,nOTU,1)))
  
  #Initialisation
  pcor<-matrix(NA, nrow = nOTU,ncol = nOTU)
  pcor.pval<-matrix(NA, nrow = nOTU,ncol = nOTU)
  colnames(pcor)<-colnames(pcor.pval)<-taxa_names
  rownames(pcor)<-rownames(pcor.pval)<-taxa_names
  
  if(!transform){
    
    if(is.null(lib_size)){
      # Creating library size by summing taxa counts per sample and time point
      lib_size <- apply(data,c(1,3),sum)
    }
    # Extract count array after excluding taxa
    count_filt<-data[,taxa_names,]
    
    #LUPINE_single and LUPINE
    if(single){
      ##**LUPINE_single with counts**##
      count_day_f<-count_filt[,,day]
      pca_res<-pca(count_day_f, ncomp=1, scale=TRUE)
      
      #loadings vector
      loading_v<- c(pca_res$loadings$X)
      #loading matrix
      loadings_m<-matrix(loading_v,ncol = 1)
      
      for(i in 1:len){
       # print(i)
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(count_day_f, center=TRUE, scale=TRUE) %*%loading_tmp
        #principal component regression on log counts+1 with an offset for library size
        r_i<-lm(log(count_day_f[,taxa1[i]]+1)~u1, offset = log(lib_size[,day]))
        r_j<-lm(log(count_day_f[,taxa2[i]]+1)~u1, offset = log(lib_size[,day]))
        #partial correlation calculation
        pcor[taxa1[i], taxa2[i]]<-pcor[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                     r_j$residuals, method = "pearson")$estimate
        pcor.pval[taxa1[i], taxa2[i]]<-pcor.pval[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                               r_j$residuals, method = "pearson")$p.value
      }
      
    }else{
      ##**LUPINE with counts**##
      count_day_f<-count_filt[,,day]
      
      ##Creating d1 to include all past data##
      day_i=day
      d1<-list(X1=count_filt[,,(day_i-1)])
      day_i=day_i-1
      index=3
      while(day_i!=1){
        d_tmp<-count_filt[,, (day_i-1)]
        d1<-list.append(d1,d_tmp)
        names(d1)[index-1]<-paste("X", index-1, sep = "")
        index=index+1
        day_i=day_i-1 
      }
        ##PLS and blockPLS##
        if(day==2){
          res_netPLS1<-pls(d1$X1, count_day_f, 
                           ncomp=1, max.iter =1000, near.zero.var = TRUE)
          
        }else{
          res_netPLS1<-block.pls(d1, count_day_f, 
                                 ncomp=1, max.iter =1000, near.zero.var = TRUE)
        }
        #loadings for nonZero variance
        loading_v<- c(res_netPLS1$loadings$Y)
        
        #loading matrix
        loadings_m<-matrix(loading_v,ncol = 1)
   
      for(i in 1:len){
        #print(i)
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-scale(count_day_f, center=TRUE, scale=TRUE) %*%loading_tmp
        #latent component regression on log counts+1 with an offset for library size
        r_i<-lm(log(count_day_f[,taxa1[i]]+1)~u1, offset = log(lib_size[,day]))
        r_j<-lm(log(count_day_f[,taxa2[i]]+1)~u1, offset = log(lib_size[,day]))
        #partial correlation calculation
        pcor[taxa1[i], taxa2[i]]<-pcor[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                     r_j$residuals, method = "pearson")$estimate
        pcor.pval[taxa1[i], taxa2[i]]<-pcor.pval[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                               r_j$residuals, method = "pearson")$p.value
      }
      
    }
    
  }else{
    
    # CLR transformation of the data
    count_full.clr <- sapply(1:dim(data)[3], function(i){
                                                clr_trans(data[,,i])}, simplify = FALSE)%>%
                      unlist( . ) %>%
                      array( ., dim = c( dim(data)[1] ,dim(data)[2] , dim(data)[3] ))
    # Extract clr array after excluding taxa
    count_filt.clr<-count_full.clr[,taxa_names,]
    
    #LUPINE_single and LUPINE
    if(single){
      ##**LUPINE_single with clr**##
      count_day_f.clr<-count_filt.clr[,,day]
      
      ##checking which taxa has zero variance to avoid errors in pca##
      zero_var_index<-which(count_day_f.clr %>% apply(.,2,var)==0)
      
      if(length(zero_var_index)>0){
        pca_res<-pca(count_day_f.clr[,-zero_var_index], ncomp=1, scale=TRUE)
        
        #loadings for nonZero variance
        loading_v<- c(pca_res$loadings$X)
        
        #Including a loading of zero to taxa with zero variance
        for(i in 1:(length(zero_var_index))) {
          if(zero_var_index[i]==1){
            loading_v <- append(0,loading_v)
          }else{
            loading_v <- append(loading_v, 0, after=(zero_var_index[i]-1))
          }
        }
  
        #loading matrix
        loadings_m<-matrix(loading_v,ncol = 1)
      } else{
        pca_res<-pca(count_day_f.clr, ncomp=1, scale=TRUE)
        
        #loadings for nonZero variance
        loading_v<- c(pca_res$loadings$X)
        
        #loading matrix
        loadings_m<-matrix(loading_v,ncol = 1)
      }
      for(i in 1:len){
        #print(i)
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-count_day_f.clr %*%loading_tmp
        #principal component regression on clr
        r_i<-lm(count_day_f.clr[,taxa1[i]]~u1)
        r_j<-lm(count_day_f.clr[,taxa2[i]]~u1)
        #partial correlation calculation
        pcor[taxa1[i], taxa2[i]]<-pcor[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                     r_j$residuals, method = "pearson")$estimate
        pcor.pval[taxa1[i], taxa2[i]]<-pcor.pval[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                               r_j$residuals, method = "pearson")$p.value
      }
    }else{
      ##**LUPINE with clr**##
      
      count_day_f.clr<-count_filt.clr[,,day]
      ##checking which taxa has zero variance to avoid errors in pls##
      zero_var_index<-which(count_day_f.clr %>% apply(.,2,var)==0)
      
      ##Creating d1 to include all past data while making sure no zero variance taxa##
      day_i=day
      zero_var_tmp<-which(count_filt.clr[,, (day_i-1)]%>% apply(.,2,var)==0)
      d1<-list(X1=if(length(zero_var_tmp)!=0)
        count_filt.clr[,-zero_var_tmp, (day_i-1)] else
          count_filt.clr[,, (day_i-1)])
      day_i=day_i-1
      index=3
      while(day_i!=1){
        zero_var_tmp<-which(count_filt.clr[,, (day_i-1)]%>% apply(.,2,var)==0)
        d_tmp<-if(length(zero_var_tmp)!=0)
          count_filt.clr[,-zero_var_tmp, (day_i-1)] else
            count_filt.clr[,, (day_i-1)]
        d1<-list.append(d1,d_tmp)
        names(d1)[index-1]<-paste("X", index-1, sep = "")
        index=index+1
        day_i=day_i-1 
      }
      
      if(length(zero_var_index)>0){
        if(day==2){
          res_netPLS1<-pls(d1$X1, count_day_f.clr[,-zero_var_index], 
                           ncomp=1, max.iter =1000)
          
        }else{
          res_netPLS1<-block.pls(d1, count_day_f.clr[,-zero_var_index], 
                                 ncomp=1, max.iter =1000)
        }
        #loadings for nonZero variance
        loading_v<- c(res_netPLS1$loadings$Y)
        
        for(i in 1:(length(zero_var_index))) {
          if(zero_var_index[i]==1){
            loading_v <- append(0,loading_v)
          }else{
            loading_v <- append(loading_v, 0, after=(zero_var_index[i]-1))
          }
        }
        
        #loading matrix
        loadings_m<-matrix(loading_v,ncol = 1)
      }else{
        if(day==2){
          res_netPLS1<-pls(d1$X1, count_day_f.clr, 
                           ncomp=1, max.iter =1000)
          
        }else{
          res_netPLS1<-block.pls(d1, count_day_f.clr, 
                                 ncomp=1, max.iter =1000)
        }
        #loadings for nonZero variance
        loading_v<- c(res_netPLS1$loadings$Y)
        
        #loading matrix
        loadings_m<-matrix(loading_v,ncol = 1)
        
      }
      for(i in 1:len){
        #print(i)
        loading_tmp<-loadings_m
        loading_tmp[c(taxa1[i], taxa2[i]),]<-0
        u1<-count_day_f.clr %*%loading_tmp
        
        r_i<-lm(count_day_f.clr[,taxa1[i]]~u1)
        r_j<-lm(count_day_f.clr[,taxa2[i]]~u1)
        pcor[taxa1[i], taxa2[i]]<-pcor[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                     r_j$residuals, method = "pearson")$estimate
        pcor.pval[taxa1[i], taxa2[i]]<-pcor.pval[taxa2[i], taxa1[i]]<-cor.test(r_i$residuals, 
                                                                               r_j$residuals, method = "pearson")$p.value
      }
    }
    
    
  }
  
  pcor.full<-matrix(NA, nrow = nOTU_total,ncol = nOTU_total)
  pcor.pval.full<-matrix(NA, nrow = nOTU_total,ncol = nOTU_total)
  colnames(pcor.full)<-colnames(pcor.pval.full)<-colnames(count_day)
  rownames(pcor.full)<-rownames(pcor.pval.full)<-colnames(count_day)
  
  common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))
  
  pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
  pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]
  
  res<-list(Estimate=pcor.full, pvalue=pcor.pval.full)
  
  return(res)
}
