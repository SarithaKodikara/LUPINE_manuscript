
library(dplyr)
library(ggplot2)
library(tidyr) 
library(grid) 
library(gridExtra) 
library(cowplot) 

load("Results/df_sen_pre_p_1.rds")
df_pca1<-df_sen_pre
load("Results/df_sen_pre_p_2.rds")
df_pca2<-df_sen_pre
load("Results/df_sen_pre_p_3.rds")
df_pca3<-df_sen_pre
###plots for PLS

df_sen_pre_all<- rbind(df_pca1,df_pca2, df_pca3)
df_sen_pre_all$Method<-factor(df_sen_pre_all$Method, 
                              levels = c( "PCA_ncomp1", "PCA_ncomp2","PCA_ncomp3"), 
                              labels = c("LUPINE_single (ncomp=1)", "LUPINE_single (ncomp=2)", "LUPINE_single (ncomp=3)"))
data_long <- gather(df_sen_pre_all, AUC, Value, ROC:PRC, factor_key=TRUE)

#df_sen_pre_fil<-filter(df_sen_pre_fil, Method %in% c("glasso","PLS_glasso","gl&PLS_gl",
# "sparcc","SpiecEasi_glasso",
#"SpiecEasi_mb","pearson"))

#df_sen_pre_fil<-df_sen_pre_fil[-c(11,12),]

pdf("Figures/n23_pca_ncomp.pdf",
    width=8, height=8)
set.seed(12345)
op <- par(mar = rep(0, 4))

data_long %>%
  ggplot(aes(x=Day, y=Value, fill=Method, color=Method)) +
  geom_boxplot() +
  ylim(0.1,0.8)+
  scale_fill_manual(values=c("deepskyblue4","turquoise4", "darkseagreen4"))+
  scale_color_manual(values=c("grey28","grey28", "grey28"))+
  facet_wrap(~AUC,labeller = labeller(AUC = 
                                        c("ROC" = "AUC-ROC",
                                          "PRC" = "AUC-PRC")
  ))+theme_bw()+
  theme(text = element_text(size=20),legend.position = 'none')

par(op)
dev.off()


#Legend

gplot <- data_long %>%
  ggplot(aes(x=Day, y=Value, fill=Method, color=Method)) +
  geom_boxplot() +
  scale_fill_manual(values=c("deepskyblue4","turquoise4", "darkseagreen4"))+
  scale_color_manual(values=c("grey28","grey28", "grey28"))+
  facet_wrap(~AUC,labeller = labeller(AUC = 
                                        c("ROC" = "AUC-ROC",
                                          "PRC" = "AUC-PRC")
  )) +theme_bw() +
  theme(legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(1, 'cm'), #change legend key height
        legend.key.width = unit(1, 'cm'), #change legend key width
        legend.title = element_text(size=14, face = "bold"), #change legend title font size
        legend.text = element_text(size=12), 
        legend.direction="horizontal") #change legend text font size

# Draw Only Legend without plot 
# Grab legend from gplot 
legend <- get_legend(gplot)                     

# Create new plot window 
grid.newpage()                               

# Draw Only legend  
grid.draw(legend)  

pdf("Figures/legend_pca_ncomp.pdf",
    width=9, height=0.5)
op <- par(mar = rep(0, 4))
grid.draw(legend) 
par(op)
dev.off()
