library(dglm)
library("cluster")
library(circlize)
library(RColorBrewer)
library(patchwork)
library(igraph)
library(ggplotify)
library(abind)
library(corpcor)
library(simstudy)
library(ggplot2)
library(gplots)
library(reshape2)
library(corpcor)
library(tidyr)
library(purrr)
library(SpiecEasi)
library(stringr)

dmlgamma <- function(x) {
  fit <- dglm(
    x~1, 
    family=Gamma(link="log"), 
    mustart=mean(x)
  )
  mu <- exp(fit$coefficients)
  shape <- exp(-fit$dispersion.fit$coefficients)
  scale <- mu/shape
  result <- c(shape, scale)
  names(result) <- c("shape", "scale")
  result
}

my_qqplot <- function(.data, .title) {
  ggplot(data = .data, mapping = aes(sample = value)) +
    stat_qq(distribution = qgamma, 
            dparams = dmlgamma(.data$value)) +
    stat_qq_line(distribution = qgamma, 
                 dparams = dmlgamma(.data$value)) +
    facet_wrap(~ name, scales = "free") +
    labs(x = "Theoretical Quantiles", y = "Sample Quantiles", title = .title)
}


load("Data_filtered/filtered_taxonomy.rds")
load("Data_filtered/OTUdata_Normal.rds")
load("Data_filtered/OTUdata_HFHS.rds")

Day1_N<-OTUdata_Normal[,,1]
Day1_H<-OTUdata_HFHS[,,1]
taxa_info_N<-filtered_taxonomy[colnames(Day1_N),]
table(taxa_info_N$X5) #1--54 o__Bacteroidales

OTUdata_Normal_bac<-OTUdata_Normal[,1:54,]
OTUdata_HFHS_bac<-OTUdata_HFHS[,1:54,]

meanOTU_n<-apply(OTUdata_Normal_bac+1, c(1,2), mean)
colnames(meanOTU_n)<-str_replace(colnames(meanOTU_n), "OTU", "Taxa")
gamma_par_n<-t(apply(meanOTU_n,2,dmlgamma))
meanOTU_n_l = pivot_longer(data.frame(meanOTU_n), cols = everything())
qqplots_n <- meanOTU_n_l %>% 
  split(.$name) %>% 
  imap(my_qqplot)


pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/qqplot1.pdf",
    width=15, height=15)
ggpubr::ggarrange(plotlist = qqplots_n, nrow = 8, ncol = 7)
dev.off()

meanOTU_h<-apply(OTUdata_HFHS_bac+1, c(1,2), mean)
colnames(meanOTU_h)<-str_replace(colnames(meanOTU_h), "OTU", "Taxa")
gamma_par_h<-t(apply(meanOTU_h,2,dmlgamma))
meanOTU_h_l = pivot_longer(data.frame(meanOTU_h), cols = everything())
qqplots_h <- meanOTU_h_l %>% 
  split(.$name) %>% 
  imap(my_qqplot)

pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/qqplot2.pdf",
    width=15, height=15)
ggpubr::ggarrange(plotlist = qqplots_h, nrow = 8, ncol = 7)
dev.off()

data_n<- apply(OTUdata_Normal, 2L, c)
spiec_n<-spiec.easi(data_n, method='mb')
binary.graph1_n <- as.matrix(spiec_n$refit$stars)[1:54, 1:54]
ig_n     <- adj2igraph(binary.graph1_n)%>%
  set_vertex_attr("name", value = paste0(1:54))

set.seed(12345)
coords <- layout_with_fr(ig_n)
plot(ig_n, vertex.size=10, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=3, layout=coords)

data_h<- apply(OTUdata_HFHS, 2L, c)
spiec_h<-spiec.easi(data_h, method='mb')
binary.graph1_h <- as.matrix(spiec_h$refit$stars)[1:54, 1:54]
ig_h     <- adj2igraph(binary.graph1_h)%>%
  set_vertex_attr("name", value = paste0(1:54))
set.seed(12345)
coords <- layout_with_fr(ig_h)
plot(ig_h, vertex.size=8,  edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=3, layout=coords, label.cex=12)

pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/simNet1_new.pdf")
set.seed(12349)
op <- par(mar = rep(0, 4))
V(ig_n)$label.cex = 1
plot(ig_n, vertex.size=8, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=2, layout=layout_with_fr(ig_n), label.cex=12)
par(op)
dev.off()


pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/simNet2_new.pdf")
set.seed(12349)
op <- par(mar = rep(0, 4))
V(ig_h)$label.cex = 1
plot(ig_h, vertex.size=8, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=2, layout=layout_with_fr(ig_h), label.cex=12)
par(op)
dev.off()




para=list(gamma_par_n=gamma_par_n,gamma_par_h=gamma_par_h,
          binary.graph1_n=binary.graph1_n, binary.graph1_h=binary.graph1_h)

saveRDS(para, file="Sim_para/para.rds")
