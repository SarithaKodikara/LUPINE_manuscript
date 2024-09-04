
library(dplyr)
library(ggplot2)
library(tidyr) 
library(grid) 
library(gridExtra) 
library(cowplot) 
library(patchwork)

load("Results/df_sen_pre_n10.rds")
df_10<-df_sen_pre%>% filter(Day==5) %>%mutate(Sample_size=10)
load("Results/df_sen_pre_n50.rds")
df_50<-df_sen_pre%>% filter(Day==5) %>%mutate(Sample_size=50)
load("Results/df_sen_pre_n100.rds")
df_100<-df_sen_pre%>% filter(Day==5) %>%mutate(Sample_size=100)
load("Results/df_sen_pre_n500.rds")
df_500<-df_sen_pre%>% filter(Day==5) %>%mutate(Sample_size=500)
load("Results/df_sen_pre_n1000.rds")
df_1000<-df_sen_pre%>% filter(Day==5) %>%mutate(Sample_size=1000)


df_sen_pre_all<- rbind(df_10,df_50, df_100,df_500,df_1000)
df_sen_pre_all$Method<-factor(df_sen_pre_all$Method, levels = c( "PCA","bPLS"), 
                              labels = c("LUPINE_single",  "LUPINE"))
data_long <- gather(df_sen_pre_all, AUC, Value, ROC:PRC, factor_key=TRUE)
data_long$Sample_size<-factor(data_long$Sample_size)

pdf("Figures/Sample_size.pdf",
    width=10, height=8)
set.seed(12345)
op <- par(mar = rep(0, 4))

data_long %>%
  ggplot(aes(x=Sample_size, y=Value, fill=Method, color=Method)) +
  geom_boxplot() +
  ylim(0,1)+
  scale_fill_manual(values=c("plum1","plum4"))+
  scale_color_manual(values=c("grey28","grey28"))+
  facet_wrap(~AUC,labeller = labeller(AUC = 
                                        c("ROC" = "AUC-ROC",
                                          "PRC" = "AUC-PRC")
  ))+theme_bw()+
  theme(text = element_text(size=20),legend.position = 'none')+
  xlab("Sample size")

par(op)
dev.off()


#Legend

gplot <- data_long %>%
  ggplot(aes(x=Day, y=Value, fill=Method, color=Method)) +
  geom_boxplot() +
  scale_fill_manual(values=c("plum1","plum4"))+
  scale_color_manual(values=c("grey28","grey28"))+
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

pdf("Figures/legend.pdf",
    width=4, height=0.5)
op <- par(mar = rep(0, 4))
grid.draw(legend) 
par(op)
dev.off()
