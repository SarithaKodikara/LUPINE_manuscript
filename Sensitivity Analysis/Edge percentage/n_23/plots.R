
library(dplyr)
library(ggplot2)
library(tidyr) 
library(grid) 
library(gridExtra) 
library(cowplot) 
library(patchwork)

load("Results/df_sen_pre_0.005.rds")
df_0.005<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.005)
load("Results/df_sen_pre_0.01.rds")
df_0.01<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.01)
load("Results/df_sen_pre_0.05.rds")
df_0.05<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.05)
load("Results/df_sen_pre_0.1.rds")
df_0.1<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.1)
load("Results/df_sen_pre_0.2.rds")
df_0.2<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.2)
load("Results/df_sen_pre_0.4.rds")
df_0.4<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.4)
load("Results/df_sen_pre_0.8.rds")
df_0.8<-df_sen_pre%>% filter(Day==5) %>%mutate(Edge_per=0.8)


df_sen_pre_all<- rbind(df_0.005,df_0.01, df_0.05,df_0.1,df_0.2,df_0.4,df_0.8)
df_sen_pre_all$Method<-factor(df_sen_pre_all$Method, levels = c( "PCA","bPLS"), 
                              labels = c("LUPINE_single",  "LUPINE"))
data_long <- gather(df_sen_pre_all, AUC, Value, ROC:PRC, factor_key=TRUE)
data_long$Edge_per<-factor(data_long$Edge_per)

pdf("Figures/Edge_per.pdf",
    width=10, height=8)
set.seed(12345)
op <- par(mar = rep(0, 4))

data_long %>%
  ggplot(aes(x=Edge_per, y=Value, fill=Method, color=Method)) +
  geom_boxplot() +
  ylim(0,1)+
  scale_fill_manual(values=c("plum1","plum4"))+
  scale_color_manual(values=c("grey28","grey28"))+
  facet_wrap(~AUC,labeller = labeller(AUC = 
                                        c("ROC" = "AUC-ROC",
                                          "PRC" = "AUC-PRC")
  ))+theme_bw()+
  theme(text = element_text(size=20),legend.position = 'none')+
  xlab("Proportion of Edges ")

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
