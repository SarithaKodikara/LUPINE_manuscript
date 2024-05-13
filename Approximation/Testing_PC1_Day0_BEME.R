library(mixOmics)
library(vegan)
library(car)
library(reshape2)
library(DescTools)
load("data_filtered/OTUdata_Normal.rds")

Y_count.day0<-OTUdata_Normal[,,1]
nOTU<-dim(Y_count.day0)[2]
len<-nOTU*(nOTU-1)/2
taxa1=unlist(lapply(1:nOTU, function(i)rep(i,(nOTU-i))))
taxa2=unlist(lapply(2:nOTU, function(i)seq(i,nOTU,1)))


pca_full<-pca(Y_count.day0, ncomp=1,scale = TRUE)
loadings_f<-matrix(c(pca_full$loadings$X),ncol = 1)


concor<-c()


for(i in 1:len){
  print(i)
  loading_tmp<-loadings_f
  loading_tmp[c(taxa1[i], taxa2[i]),]<-0
  u1_full<-scale(Y_count.day0) %*%loading_tmp
  
  pca_apro<-pca(Y_count.day0[,-c(taxa1[i], taxa2[i])], ncomp=1, scale = TRUE)
  u1_apro<-pca_apro$variates$X
  
  concor[i] <-CCC(u1_full[,1],u1_apro[,1])$rho.c[[1]]
}

df_concor<-data.frame(concor)

pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/supp_1.pdf",
    width=8, height=8)
set.seed(12345)
op <- par(mar = rep(0, 4))
df_concor %>%
  ggplot(aes( y=(concor))) +
  geom_boxplot(color="blue", fill="skyblue", alpha=0.2) +
  theme(text = element_text(size=30),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())+ 
  ylab("Concordance correlation")
par(op)
dev.off()
