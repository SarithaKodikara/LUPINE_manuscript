
nworkers <- parallel::detectCores() ## number of cpus to use in paralle processing
cat("\nUsing ", nworkers, " CPUs ...\n")
BPPARAM <- BiocParallel::MulticoreParam(workers = nworkers)

source("LUPINE_function.R")
load("data_modified/data_MDE.rds")
load("data_modified/data_PPT.rds")

# outer lapply for iter
inf_upper<-BiocParallel::bplapply(
  2,
  function(day){
    # second inner lapply
    res_2<-lapply(
      c('PPT','MDE'),
      function(m){
        if(m=='MDE'){
            res<-LUPINE(data=data_MDE, day=day,  transformed=TRUE, single = FALSE)
            save(res, file=paste0("Results/",m,"/LUPINE_Day",day,".rdata"))
          
        }else{
          res<-LUPINE(data=data_PPT, day=day,  transformed=TRUE, single = FALSE)
          save(res, file=paste0("Results/",m,"/LUPINE_Day",day,".rdata"))
        }
      })
  })