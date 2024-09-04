library(parallel)
nworkers <- parallel::detectCores() ## number of cpus to use in paralle processing

source("auc_roc_prc.R")

runFunction<-function(iter){
    res_1<-lapply(
      2:10,
      function(day){
        # second inner lapply
        res_2<-lapply(
          c('bPLS_clr'), 
          function(m){
            auc_res( iter, day, method=m)})
        # concatenate second inner rows
        do.call(rbind,res_2)
        
      })
    # concatenate first inner rows
    return(do.call(rbind,res_1))
}

# outer lapply for iter
inf_upper<-mclapply(1:50, runFunction, mc.cores = nworkers-2)


df_sen_pre<-do.call(rbind, inf_upper)
df_sen_pre$Day<-factor(df_sen_pre$Day)
save(df_sen_pre,file = "Results/df_sen_pre_b_clr.rds")
