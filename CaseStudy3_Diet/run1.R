
nworkers <- parallel::detectCores() ## number of cpus to use in paralle processing
cat("\nUsing ", nworkers, " CPUs ...\n")
BPPARAM <- BiocParallel::MulticoreParam(workers = nworkers)

source("LUPINE_function.R")
load("data_modified/OTUdata_plant.rds")
load("data_modified/low_plant.rds")
load("data_modified/OTUdata_animal.rds")
load("data_modified/low_animal.rds")

# outer lapply for iter
inf_upper<-BiocParallel::bplapply(
  2:15,
  function(day){
    # second inner lapply
    res_2<-lapply(
      c('Plant',"Animal"),
      function(m){
        if(m=='Plant'){
            res<-LUPINE(data=OTUdata_plant, day=day, excluded_taxa=low_plant, transformed=TRUE, single = FALSE)
            save(res, file=paste0("Results/",m,"/LUPINE_Day",day,".rdata"))
          
        }else{
          res<-LUPINE(data=OTUdata_animal, day=day, excluded_taxa=low_animal, transformed=TRUE, single = FALSE)
          save(res, file=paste0("Results/",m,"/LUPINE_Day",day,".rdata"))
        }
      })
  })