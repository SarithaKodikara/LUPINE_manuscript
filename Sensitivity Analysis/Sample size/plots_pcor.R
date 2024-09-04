
load("Sim_Data/pcor_e0.005.rds")
pcor_0.005<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.005<- data.frame(pcor=pcor_0.005[pcor_0.005!=0], Edge_per=0.005)
load("Sim_Data/pcor_e0.01.rds")
pcor_0.01<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.01<- data.frame(pcor=pcor_0.01[pcor_0.01!=0], Edge_per=0.01)
load("Sim_Data/pcor_e0.05.rds")
pcor_0.05<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.05<- data.frame(pcor=pcor_0.05[pcor_0.05!=0], Edge_per=0.05)
load("Sim_Data/pcor_e0.1.rds")
pcor_0.1<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.1<- data.frame(pcor=pcor_0.1[pcor_0.1!=0], Edge_per=0.1)
load("Sim_Data/pcor_e0.2.rds")
pcor_0.2<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.2<- data.frame(pcor=pcor_0.2[pcor_0.2!=0], Edge_per=0.2)
load("Sim_Data/pcor_e0.4.rds")
pcor_0.4<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.4<- data.frame(pcor=pcor_0.4[pcor_0.4!=0], Edge_per=0.4)
load("Sim_Data/pcor_e0.8.rds")
pcor_0.8<-abs(unlist(lapply(pcorData, function(mat) mat[upper.tri(mat)])))
df_0.8<- data.frame(pcor=pcor_0.8[pcor_0.8!=0], Edge_per=0.8)

df<- rbind(df_0.005,df_0.01, df_0.05,df_0.1,df_0.2,df_0.4,df_0.8)
df$Edge_per<-factor(df$Edge_per)

pdf("Figures/pcor.pdf",
    width=6, height=8)
set.seed(12345)
op <- par(mar = rep(0, 4))

df %>%
  ggplot(aes(x=Edge_per, y=pcor, fill=Edge_per, color=Edge_per)) +
  geom_boxplot() +
  ylim(0,1)+
  scale_fill_manual(values=c("cornsilk","antiquewhite1","peachpuff","bisque2", "peachpuff3", "lightsalmon4", "peachpuff4"))+
  scale_color_manual(values=c("grey28","grey28","grey28","grey28","grey28","grey28","grey28"))+theme_bw()+
  theme(text = element_text(size=20),legend.position = 'none')+
  xlab("Percentage of Edges ")+
  ylab("Absolute non-zero off-diagonal partial correlations")

par(op)
dev.off()

