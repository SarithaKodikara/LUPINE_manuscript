
nworkers <- parallel::detectCores() ## number of cpus to use in paralle processing
cat("\nUsing ", nworkers, " CPUs ...\n")
BPPARAM <- BiocParallel::MulticoreParam(workers = nworkers)

source("LUPINE_function.R")

load("data_filtered/OTUdata_Normal.rds")
load("data_filtered/Lib_Normal.rds")
load("data_filtered/low_Normal.rds")
# outer lapply for iter
inf_upper<-BiocParallel::bplapply(
      1:4,
      function(day){
        # second inner lapply
        res_2<-lapply(
          c('LUPINE', 'LUPINE_single'), 
          function(m){
            if(m=='LUPINE'){
              if(day>1){
                res<-LUPINE(data=OTUdata_Normal, day=day, excluded_taxa=low_Normal, lib_size = Lib_Normal, single = FALSE)
                save(res, file=paste0("Results/Normal/",m,"_Day",day,".rdata"))
              }
              
            }else{
              res<-LUPINE(data=OTUdata_Normal, day=day, excluded_taxa=low_Normal, lib_size = Lib_Normal, single = TRUE)
              save(res, file=paste0("Results/Normal/",m,"_Day",day,".rdata"))
            }
              })
      })


load("data_filtered/OTUdata_HFHS.rds")
load("data_filtered/Lib_HFHS.rds")
load("data_filtered/low_HFHS.rds")
# outer lapply for iter
inf_upper<-BiocParallel::bplapply(
  1:4,
  function(day){
    # second inner lapply
    res_2<-lapply(
      c('LUPINE', 'LUPINE_single'), 
      function(m){
        if(m=='LUPINE'){
          if(day>1){
            res<-LUPINE(data=OTUdata_HFHS, day=day, excluded_taxa=low_HFHS, lib_size = Lib_HFHS, single = FALSE)
            save(res, file=paste0("Results/HFHS/",m,"_Day",day,".rdata"))
          }
          
        }else{
          res<-LUPINE(data=OTUdata_HFHS, day=day, excluded_taxa=low_HFHS, lib_size = Lib_HFHS, single = TRUE)
          save(res, file=paste0("Results/HFHS/",m,"_Day",day,".rdata"))
        }
      })
  })
