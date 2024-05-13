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


load("Data_filtered/filtered_taxonomy.rds")
load("Data_filtered/OTUdata_Normal.rds")
load("Data_filtered/OTUdata_HFHS.rds")

Day1_N<-OTUdata_Normal[,,1]
Day1_H<-OTUdata_HFHS[,,1]
taxa_info_N<-filtered_taxonomy[colnames(Day1_N),]
table(taxa_info_N$X5) #1--54 o__Bacteroidales

#Calculate Bacteroidales mean abundance on day 1
meanOTU_n<-colMeans(Day1_N[,1:54])

spiec_n<-spiec.easi(Day1_N, method = "mb")
binary.graph1_n <- as.matrix(spiec_n$refit$stars)[1:54, 1:54]
ig_n     <- adj2igraph(binary.graph1_n)%>%
  set_vertex_attr("name", value = paste0(1:54))

spiec_h<-spiec.easi(Day1_H, method = "mb")
binary.graph1_h <- as.matrix(spiec_h$refit$stars)[1:54, 1:54]
ig_h     <- adj2igraph(binary.graph1_h)%>%
  set_vertex_attr("name", value = paste0(1:54))

set.seed(0178)
coords <- layout_with_dh(ig_n)
plot(ig_n, vertex.size=8, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=3, layout=coords)


set.seed(0178)
coords <- layout_nicely(ig_h)
plot(ig_h, vertex.size=10, vertex.size=8, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=3, layout=coords)

pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/simNet1_new.pdf")
set.seed(0178)
op <- par(mar = rep(0, 4))
V(ig_n)$label.cex = 1.5
plot(ig_n, vertex.size=8, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=3, layout=layout_nicely(ig_h))
par(op)
dev.off()


pdf("/Users/skodikara/Library/CloudStorage/OneDrive-TheUniversityofMelbourne/Documents/ALL/PostDoc/TempoNet/simNet2_new.pdf")
set.seed(0178)
op <- par(mar = rep(0, 4))
V(ig_h)$label.cex = 1.5
plot(ig_h, vertex.size=10, edge.color="lightsteelblue3", vertex.color="tomato",
     edge.width=3, layout=layout_nicely(ig_h), label.cex=15)
par(op)
dev.off()




para=list(mean_otu=meanOTU_n,
          binary.graph1_n=binary.graph1_n, binary.graph1_h=binary.graph1_h)

saveRDS(para, file="Sim_para/para.rds")
