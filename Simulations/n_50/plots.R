
library(dplyr)
library(ggplot2)
library(tidyr)
load("Results/df_sen_pre_m.rds")
df_mb<-df_sen_pre
load("Results/df_sen_pre_g.rds")
df_gl<-df_sen_pre
load("Results/df_sen_pre_s.rds")
df_sp<-df_sen_pre
load("Results/df_sen_pre_p.rds")
df_pca<-df_sen_pre
load("Results/df_sen_pre_b.rds")
df_bpls<-df_sen_pre

###plots for PLS

df_sen_pre_all<- rbind(df_mb,df_gl,df_sp,
                       df_pca, df_bpls)
df_sen_pre_all$Method<-factor(df_sen_pre_all$Method, 
                              levels = c("sparcc","SpiecEasi_glasso",
                              "SpiecEasi_mb", "PCA", "bPLS"), 
                              labels = c("sparcc","SpiecEasi_glasso",
                                         "SpiecEasi_mb", "TempoNet_single", "TempoNet_longi"))
data_long <- gather(df_sen_pre_all, AUC, Value, ROC:PRC, factor_key=TRUE)

#df_sen_pre_fil<-filter(df_sen_pre_fil, Method %in% c("glasso","PLS_glasso","gl&PLS_gl",
# "sparcc","SpiecEasi_glasso",
#"SpiecEasi_mb","pearson"))

#df_sen_pre_fil<-df_sen_pre_fil[-c(11,12),]

pdf("Figures/sim_n50_n.pdf",
    width=36, height=26)
set.seed(12345)
op <- par(mar = rep(0, 4))

data_long %>%
  ggplot(aes(x=Day, y=Value, fill=Method, color=Method)) +
  geom_boxplot() +
  ylim(0,1)+
  scale_fill_manual(values=c("grey", "#CCFF33", "#33FF99","plum1","magenta"))+
  scale_color_manual(values=c("grey28","olivedrab", "green4","plum4", "magenta4"))+
  facet_wrap(~AUC,labeller = labeller(AUC = 
                                       c("ROC" = "AUC-ROC",
                                         "PRC" = "AUC-PRC")
  )) +theme_bw()+
  theme(text = element_text(size=20)) + theme(legend.position = 'none')
par(op)
dev.off()
