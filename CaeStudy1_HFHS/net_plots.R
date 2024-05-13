library(RColorBrewer)
library(patchwork)
library(igraph)
library(dplyr)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(SpiecEasi)
library(ggraph)   
library(graphlayouts)
library(tidyverse)
library(tidygraph)

load("data_filtered/filtered_taxonomy.rds")
load("data_filtered/OTUdata_Normal.rds")
load("data_filtered/OTUdata_HFHS.rds")


netPlotHFHS<-function(data, taxonomy, diet, day, method){
  Day0<-data[,,1]
  Colors=brewer.pal(9,"Greys")
  Colors=colorRampPalette(Colors)(100)
  
  col = colorRamp2(seq(0,1,length=100), colorRampPalette(Colors)(100))
  
  
  taxa_info<-taxonomy[colnames(Day0),]
  taxa_info$X5<-factor(taxa_info$X5)
  col2 = list(Order = c(" o__Bacteroidales" = "green", " o__Bifidobacteriales" = "gray", 
                        " o__Burkholderiales"= "darkgreen", " o__CW040"="pink", " o__Deferribacterales"= "tomato",
                        " o__Clostridiales" = "darkred"," o__Coriobacteriales" = "firebrick2", 
                        " o__Desulfovibrionales" = "orange"," o__Erysipelotrichales"="blue", 
                        " o__Lactobacillales"="purple", " o__YS2"="hotpink", " o__Verrucomicrobiales"="lightblue"
  ))
  # Create the heatmap annotation
  ha <- HeatmapAnnotation(
    Order = taxa_info$X5,
    col = col2, show_legend = FALSE
  )
  
  
  row_ha = rowAnnotation( Order = taxa_info$X5, col = col2)
  
  load(paste0("results/",diet,"/",method,"_Day", day, ".rdata"))
  net<-(res$pvalue<0.05)*1
  net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
  
  
  Heatmap(net, cluster_rows = F, cluster_columns = F, col = col,
          top_annotation = ha,left_annotation = row_ha,
          row_split =rep(LETTERS[1:12],c(summary(taxa_info$X5))),
          column_split =rep(LETTERS[1:12],c(summary(taxa_info$X5))),
          column_names_gp= gpar(fontsize = 0), row_names_gp= gpar(fontsize = 0),
          row_title =NULL, column_title = NULL)
  
  
  g <- graph.adjacency(net, mode="undirected", weighted=NULL)
  
  # provide some names
  V(g)$name <- 1:vcount(g)
  
  
  # plot using ggraph
  graph_tbl <- g %>% 
    as_tbl_graph() %>% 
    activate(nodes) %>% 
    mutate(degree  = centrality_degree()) %>%
    mutate(community = as.factor(c(rep("green",54), rep("gray",2), rep("darkgreen",2),
                                   rep("darkred",115), rep("firebrick2",7), rep("pink",1),
                                   rep("tomato",1), rep("orange",3), rep("blue",15),
                                   rep("purple",9), rep("hotpink",1), rep("lightblue",2)))) 
  
  
  layout <- create_layout(graph_tbl, layout = 'igraph', algorithm = 'sphere')
  layout$x[1:54]<-layout$x[1:54]
  
  layout$x[55:56]<-layout$x[55:56]+0.8
  layout$y[55:56]<-layout$y[55:56]+3.5
  
  layout$x[57:58]<-layout$x[57:58]+0.5
  layout$y[57:58]<-layout$y[57:58]+3
  
  layout$y[59:173]<-layout$y[59:173]+2.5
  
  layout$x[174:180]<-layout$x[174:180]+2
  
  layout$y[181]<-layout$y[181]+2.5
  layout$x[181]<-layout$x[181]+2.5
  
  layout$y[182]<-layout$y[182]+2
  layout$x[182]<-layout$x[182]-1
  
  layout$y[183:185]<-layout$y[183:185]+2
  layout$x[183:185]<-layout$x[183:185]-1
  
  layout$y[186:200]<-layout$y[186:200]+2
  layout$x[186:200]<-layout$x[186:200]+2
  
  layout$y[201:209]<-layout$y[201:209]-0.5
  layout$x[201:209]<-layout$x[201:209]+1.6
  
  layout$x[210]<-layout$x[210]-1.5
  
  layout$x[211:212]<-layout$x[211:212]-1
  layout$y[211:212]<-layout$y[211:212]+1
  
  p<-ggraph(layout) +
    geom_edge_fan(
      aes(color = as.factor(from), alpha = 0.2),
      show.legend = F
    ) +theme_graph(background = "white")+
    geom_node_point(
      aes(size = degree, color = as.factor(name)),
      show.legend = F
    ) +
    scale_color_manual(
      limits = as.factor(layout$name),
      values = c(rep("green",54), rep("gray",2), rep("darkgreen",2),
                 rep("darkred",115), rep("firebrick2",7), rep("pink",1),
                 rep("tomato",1), rep("orange",3), rep("blue",15),
                 rep("purple",9), rep("hotpink",1), rep("lightblue",2))
    ) +
    scale_edge_color_manual(
      limits = as.factor(layout$name),
      values = c(rep("green",54), rep("gray",2), rep("darkgreen",2),
                 rep("darkred",115), rep("firebrick2",7), rep("pink",1),
                 rep("tomato",1), rep("orange",3), rep("blue",15),
                 rep("purple",9), rep("hotpink",1), rep("lightblue",2))
    ) 
  
  return(p)
  
}

netPlotHFHS(OTUdata_Normal, filtered_taxonomy, "Normal", 2, "LUPINE")
netPlotHFHS(OTUdata_Normal, filtered_taxonomy, "Normal", 3, "LUPINE")
netPlotHFHS(OTUdata_Normal, filtered_taxonomy, "Normal", 4, "LUPINE")

netPlotHFHS(OTUdata_HFHS, filtered_taxonomy, "HFHS", 2, "LUPINE")
netPlotHFHS(OTUdata_HFHS, filtered_taxonomy, "HFHS", 3, "LUPINE")
netPlotHFHS(OTUdata_HFHS, filtered_taxonomy, "HFHS", 4, "LUPINE")
  
  
#Legend
pdf("Figures/BEME_legend_net.pdf",
    width=9.5, height=1.6)
op <- par(mar = rep(0, 4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
MyOrder = matrix(1:12, nrow = 4, ncol = 3, byrow = T)
legend("topleft", legend =c(" o__Bacteroidales" , " o__Bifidobacteriales" , 
                            " o__Burkholderiales", " o__Clostridiales" ,
                            " o__Coriobacteriales" ,
                            " o__CW040", " o__Deferribacterales",
                            " o__Desulfovibrionales" ," o__Erysipelotrichales", 
                            " o__Lactobacillales", " o__Verrucomicrobiales", " o__YS2")[MyOrder],
       pch=16, pt.cex=2, cex=1.5, 
       col = c("green","gray","darkgreen",
               "darkred","firebrick2","pink",
               "tomato","orange","blue",
               "purple","hotpink","lightblue")[MyOrder],
       box.lwd = 0.1,box.col = "white",bg = "white",ncol=3)
par(op)
dev.off()
