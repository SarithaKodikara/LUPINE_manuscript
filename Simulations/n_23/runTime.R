library(dplyr)
library(scales)

runTime<-function(Iter,  Day, Method="bPLS"){
    load(paste0("Results/",Method,"/iter",Iter,"_d",Day,".rdata"))
  
      runT<-res$RunTime[3]/60
  
    return(runT)
  
}


df_runT <- tibble(Iter = rep(rep(1:50,each=9),5), Day= rep(2:10,250), Method = rep(c("SpiecEasi_mb","sparcc", 
                                                                                     "SpiecEasi_glasso",
                                                                                     "PCA", "bPLS"),each=450))
df_runT %<>% rowwise%>% 
          mutate(elapsedTime= runTime(Iter, Day, Method))%>%
          data.frame()
df_runT$Method<-factor(df_runT$Method, 
                              levels = c("sparcc","SpiecEasi_glasso",
                                         "SpiecEasi_mb", "PCA", "bPLS"), 
                              labels = c("sparcc","SpiecEasi_glasso",
                                         "SpiecEasi_mb", "LUPINE_single", "LUPINE_longi"))

df_runT$Day<-factor(df_runT$Day)
df_runT$Dummy<-rep("Elapsed time (minutes)",2250)

pdf("Figures/sim_n23_runTime_n.pdf",
    width=20, height=16)
set.seed(12345)
op <- par(mar = rep(0, 4))

df_runT %>%
  ggplot(aes(x=Day, y=(elapsedTime), fill=Method, color=Method)) +
  geom_boxplot() + 
  scale_fill_manual(values=c("grey", "#CCFF33", "#33FF99","plum1","magenta"))+
  scale_color_manual(values=c("grey28","olivedrab", "green4","plum4", "magenta4"))+
  facet_wrap(~Dummy) +
  ylab("")+theme_bw()+ 
  theme(text = element_text(size=18),legend.position = 'none')
par(op)
dev.off()
