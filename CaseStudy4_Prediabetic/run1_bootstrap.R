
nworkers <- parallel::detectCores() ## number of cpus to use in paralle processing
cat("\nUsing ", nworkers, " CPUs ...\n")
BPPARAM <- BiocParallel::MulticoreParam(workers = nworkers)

source("LUPINE_function.R")
load("data_modified/data_MDE.rds")
load("data_modified/data_PPT.rds")
set.seed(1234)

# outer lapply for iter
inf_upper<-BiocParallel::bplapply(
  c('PPT','MDE'),
  function(m){
    # second inner lapply
    res_2<-lapply(
      1:1000,
      function(iter){
        cat(paste("\n","iter=", iter, "\n"), file="iteration.txt", append = TRUE)
        
        day=2
        if(m=='MDE'){
          new_MDE<- data_MDE[sample(1:dim(data_MDE)[1], replace = T),,]
          rownames(new_MDE)<-rownames(data_MDE)
          net<-LUPINE(data=new_MDE, day=day,  transformed=TRUE, single = FALSE)$pvalue
          net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
        }else{
          new_PPT <- data_PPT[sample(1:dim(data_PPT)[1], replace = T),,]
          rownames(new_PPT)<-rownames(data_PPT)
          net<-LUPINE(data=new_PPT, day=day,  transformed=TRUE, single = FALSE)$pvalue
          net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
        }
      })
  })

save(inf_upper, file="Results/LUPINE_boot.rdata")
