
nworkers <- parallel::detectCores() ## number of cpus to use in paralle processing
cat("\nUsing ", nworkers, " CPUs ...\n")
BPPARAM <- BiocParallel::MulticoreParam(workers = nworkers)

source("LUPINE_function.R")
load("data_filtered/OTUdata.rds")
load("data_filtered/Lib_size.rds")
load("data_filtered/OTU_l.abundance.rds")

# outer lapply for iter
inf_upper<-BiocParallel::bplapply(
      1:10,
      function(day){
        # second inner lapply
        res_2<-lapply(
          c('LUPINE'),
          function(m){
            if(m=='LUPINE'){
              if(day>1){
                res<-LUPINE(data=OTUdata_array, day=day, excluded_taxa=OTU_l.abundance, lib_size = Lib_size, single = FALSE)
                save(res, file=paste0("Results/",m,"_Day",day,".rdata"))
              }
              
            }else{
              res<-LUPINE(data=OTUdata_array, day=day, excluded_taxa=OTU_l.abundance, lib_size = Lib_size, single = TRUE)
              save(res, file=paste0("Results/",m,"_Day",day,".rdata"))
            }
          })
      })



