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
load("data_filtered/taxonomy.rds")
load("data_filtered/OTUdata.rds")

netPlotVRE<-function(data, taxonomy, day, method){
  load(paste0("Results/", method,"_Day", day,".rdata"))
  net<-(res$pvalue<0.05)*1
  net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
  g <- graph.adjacency(net, mode="undirected", weighted=NULL)
  # provide some names
  V(g)$name <- 1:vcount(g)
  
  
  Day0<-OTUdata_array[,,1]
  taxa_info<-taxanomy_filter_ordered[colnames(Day0),]
  taxa_info$V6<-factor(taxa_info$V6)
  col2 = list(Order = c(" o__Anaeroplasmatales"="grey48",
                        " o__Bacillales"="pink",
                        " o__Bacteroidales" = "green",
                        " o__Clostridiales" = "darkred",
                        " o__Deferribacterales"= "orange",
                        " o__Enterobacteriales"="red",
                        " o__Erysipelotrichales"="deepskyblue",
                        " o__Lactobacillales"="purple",
                        " o__RF39"="hotpink",
                        " o__Rickettsiales"= "darkgreen",
                        " o__Streptophyta"="yellow",
                        " o__Turicibacterales"="tomato",
                        " o__Verrucomicrobiales"="blue",
                        " o__"="black"))
  
  row_ha = rowAnnotation( Order = taxa_info$V6, col = col2)
  
  # Create the heatmap annotation
  ha <- HeatmapAnnotation(
    Order = taxa_info$V6,
    col = col2, show_legend = FALSE
  )
  
  # # plot using ggraph
  graph_tbl <- g %>%
    as_tbl_graph() %>%
    activate(nodes) %>%
    mutate(degree  = centrality_degree()) %>%
    mutate(community = as.factor(rep(c("black","grey48","pink", "green", "darkred", "orange","red","deepskyblue",
                                       "purple", "hotpink","darkgreen","yellow","tomato", "blue"),
                                     c(summary(factor(taxanomy_filter_ordered$V6))))))
  
  # Colors=brewer.pal(9,"Greys")
  # Colors=colorRampPalette(Colors)(100)
  # 
  # col = colorRamp2(seq(0,1,length=100), colorRampPalette(Colors)(100))
  # Heatmap(net, cluster_rows = F, cluster_columns = F, col = col,
  #         top_annotation = ha,left_annotation = row_ha)
  
  layout <- create_layout(graph_tbl, layout = 'igraph', algorithm = 'sphere')
  layout$x[1:14]<-layout$x[1:14]-0.5
  layout$y[1:14]<-layout$y[1:14]-3
  
  layout$x[15:16]<-layout$x[15:16]+1.2
  layout$y[15:16]<-layout$y[15:16]+1.1
  
  layout$x[17:19]<-layout$x[17:19]+0.2
  layout$y[17:19]<-layout$y[17:19]+1.1
  
  layout$x[20:27]<-layout$x[20:27]-0.5
  layout$y[20:27]<-layout$y[20:27]+1
  
  layout$x[28:111]<-layout$x[28:111]+0.5
  layout$y[28:111]<-layout$y[28:111]-1
  
  layout$x[112]<-layout$x[112]-0.8
  layout$y[112]<-layout$y[112]-1.5
  
  layout$x[113]<-layout$x[113]-0.8
  layout$y[113]<-layout$y[113]-1.5
  
  layout$x[114]<-layout$x[114]-0.8
  layout$y[114]<-layout$y[114]-1.5
  
  layout$x[115:119]<-layout$x[115:119]+1.5
  layout$y[115:119]<-layout$y[115:119]-1.5
  
  layout$x[120:122]<-layout$x[120:122]+2.1
  layout$y[120:122]<-layout$y[120:122]-0.5
  
  layout$x[123]<-layout$x[123]-0.8
  layout$y[123]<-layout$y[123]-1.5
  
  layout$x[124]<-layout$x[124]-0.8
  layout$y[124]<-layout$y[124]-1.5
  
  layout$x[125]<-layout$x[125]-0.8
  layout$y[125]<-layout$y[125]-1.5
  
  layout$x[126]<-layout$x[126]-0.8
  layout$y[126]<-layout$y[126]-1.5
  
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
      values = rep(c("black","grey48","pink", "green", "darkred", "orange","red","deepskyblue", 
                     "purple", "hotpink","darkgreen","yellow","tomato", "blue"),
                   c(summary(factor(taxanomy_filter_ordered$V6))))
    ) +
    scale_edge_color_manual(
      limits = as.factor(layout$name),
      values = rep(c("black","grey48","pink", "green", "darkred", "orange","red","deepskyblue", 
                     "purple", "hotpink","darkgreen","yellow","tomato", "blue"),
                   c(summary(factor(taxanomy_filter_ordered$V6))))
    ) 
  
  return(p)
  
}

plots<-sapply(2:10, function(i){netPlotVRE(OTUdata_array, taxanomy_filter_ordered,  i, "LUPINE")}, simplify = FALSE)
plots[[1]]
#.....
plots[[9]]

#Legend
pdf("Figures/VRE_legend_net.pdf",
    width=9.5, height=1.7)
op <- par(mar = rep(0, 4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
MyOrder = matrix(1:15, nrow = 5, ncol = 3, byrow = T)
legend("topleft", legend =c(" o__Anaeroplasmatales"," o__Bacillales",
                            " o__Bacteroidales" ," o__Clostridiales" ,
                            " o__Deferribacterales"," o__Enterobacteriales",
                            " o__Erysipelotrichales"," o__Lactobacillales",
                            " o__RF39"," o__Rickettsiales"," o__Streptophyta",
                            " o__Turicibacterales"," o__Verrucomicrobiales",
                            " o__")[MyOrder],
       pch=16, pt.cex=2, cex=1.5, 
       col = c("grey48","pink", "green", "darkred",
               "orange","red","deepskyblue","purple",
               "hotpink", "darkgreen","yellow","tomato",
               "blue","black")[MyOrder],
       box.lwd = 0.1,box.col = "white",bg = "white",ncol=3)
par(op)
dev.off()




