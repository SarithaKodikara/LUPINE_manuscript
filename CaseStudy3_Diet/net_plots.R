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


load("Results/Animal/LUPINE_Day15.rdata")
load("data_modified/taxonomy_fil.rds")

netPlotDiet<-function(data, taxonomy, diet,){
  
}



set.seed(123456)
#4, 9, 15

net<-(res$pvalue<0.05)*1
net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})

g <- graph_from_adjacency_matrix(net, mode="undirected", weighted=NULL)

# provide some names
V(g)$name <- 1:vcount(g)



# plot using ggraph
graph_tbl <- g %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(degree  = centrality_degree()) %>%
  mutate(community = as.factor(rep(
    c( "black","yellow", "green", "gray", 
       "darkgreen",  "pink", "hotpink",  "royalblue", 
       "darkred", "firebrick2",  "bisque4", "orange" , 
       "cyan3", "blue",  "darkolivegreen3", "plum1", 
       "coral2", "purple",  "aquamarine3", "khaki4", 
       "indianred", "mediumvioletred",  "tan4",
       "turquoise4",  "wheat", "springgreen", 
       "steelblue4", "lightblue", "goldenrod"),
    c(summary(factor(taxanomy_filter_ordered$Order))))))


##############
layout1<-as.matrix(create_layout(induced_subgraph(graph_tbl, 1), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout2<-as.matrix(create_layout(induced_subgraph(graph_tbl, 2:3), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2) 
layout3<-as.matrix(create_layout(induced_subgraph(graph_tbl,4:124), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout4<-as.matrix(create_layout(induced_subgraph(graph_tbl,125:126), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout5<-as.matrix(create_layout(induced_subgraph(graph_tbl,127:133), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout6<-as.matrix(create_layout(induced_subgraph(graph_tbl,134:136), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout7<-as.matrix(create_layout(induced_subgraph(graph_tbl,137:146), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout8<-as.matrix(create_layout(induced_subgraph(graph_tbl,147), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout9<-as.matrix(create_layout(induced_subgraph(graph_tbl,148:151), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout10<-as.matrix(create_layout(induced_subgraph(graph_tbl,152:153), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout11<-as.matrix(create_layout(induced_subgraph(graph_tbl,154), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout12<-as.matrix(create_layout(induced_subgraph(graph_tbl,155:156), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout13<-as.matrix(create_layout(induced_subgraph(graph_tbl,157:164), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout14<-as.matrix(create_layout(induced_subgraph(graph_tbl,165:172), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout15<-as.matrix(create_layout(induced_subgraph(graph_tbl,173:174), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout16<-as.matrix(create_layout(induced_subgraph(graph_tbl,175:176), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout17<-as.matrix(create_layout(induced_subgraph(graph_tbl,177:236), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout18<-as.matrix(create_layout(induced_subgraph(graph_tbl,237:239), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout19<-as.matrix(create_layout(induced_subgraph(graph_tbl,240), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout20<-as.matrix(create_layout(induced_subgraph(graph_tbl,241), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout21<-as.matrix(create_layout(induced_subgraph(graph_tbl,242:289), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout22<-as.matrix(create_layout(induced_subgraph(graph_tbl,290:293), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout23<-as.matrix(create_layout(induced_subgraph(graph_tbl,294), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout24<-as.matrix(create_layout(induced_subgraph(graph_tbl,295:298), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout25<-as.matrix(create_layout(induced_subgraph(graph_tbl,299), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout26<-as.matrix(create_layout(induced_subgraph(graph_tbl,300:303), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout27<-as.matrix(create_layout(induced_subgraph(graph_tbl,304:308), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout28<-as.matrix(create_layout(induced_subgraph(graph_tbl,309:310), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout29<-as.matrix(create_layout(induced_subgraph(graph_tbl,311), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)

layout0 <- list(layout1, layout2, layout3, layout4,
              layout5,layout6,layout7,layout8,
              layout9,layout10,layout11,layout12,
              layout13,layout14,layout15,layout16,
              layout17,layout18,layout19,layout20,
              layout21,layout22,layout23,layout24,
              layout25,layout26,layout27,layout28,
              layout29)
graphs<-list(induced_subgraph(graph_tbl, 1),
             induced_subgraph(graph_tbl, 2:3),
             induced_subgraph(graph_tbl, 4:124),
             induced_subgraph(graph_tbl, 125:126),
             induced_subgraph(graph_tbl, 127:133),
             induced_subgraph(graph_tbl, 134:136),
             induced_subgraph(graph_tbl, 137:146),
             induced_subgraph(graph_tbl,147),induced_subgraph(graph_tbl,148:151),
             induced_subgraph(graph_tbl,152:153),induced_subgraph(graph_tbl,154),
             induced_subgraph(graph_tbl,155:156),induced_subgraph(graph_tbl,157:164),
             induced_subgraph(graph_tbl,165:172),induced_subgraph(graph_tbl,173:174),
             induced_subgraph(graph_tbl,175:176),induced_subgraph(graph_tbl,177:236),
             induced_subgraph(graph_tbl,237:239),induced_subgraph(graph_tbl,240),
             induced_subgraph(graph_tbl,241),induced_subgraph(graph_tbl,242:289),
             induced_subgraph(graph_tbl,290:293),induced_subgraph(graph_tbl,294),
             induced_subgraph(graph_tbl,295:298),induced_subgraph(graph_tbl,299),
             induced_subgraph(graph_tbl,300:303),induced_subgraph(graph_tbl,304:308),
             induced_subgraph(graph_tbl,309:310),induced_subgraph(graph_tbl,311))
lay <- merge_coords(graphs, layout0)

layout <- create_layout(graph_tbl, layout = 'igraph', algorithm = 'circle')
layout$x<-lay[,1]+layout$x*2
layout$y<-lay[,2]+layout$y*2

op <- par(mar = rep(0, 4))
g1<-ggraph(layout) +
  geom_edge_fan(
    aes(color = as.factor(from)),
    show.legend = F, alpha = 0.5,width=0.5
  ) +theme_graph(background = "white")+
  geom_node_point(
    aes(fill = as.factor(name), color=as.factor(name),
        size=degree),shape=21,
    show.legend = F
  ) + scale_color_manual(
    values = rep(alpha("black",0.5),311)
  ) +
  scale_edge_color_manual(
    limits = as.factor(layout$name),
    values = rep(c( "black","yellow", "green", "gray", 
                    "darkgreen",  "pink", "hotpink",  "royalblue", 
                    "darkred", "firebrick2",  "bisque4", "orange" , 
                    "cyan3", "blue",  "darkolivegreen3", "plum1", 
                    "coral2", "purple",  "aquamarine3", "khaki4", 
                    "indianred", "mediumvioletred",  "tan4",
                    "turquoise4",  "wheat", "springgreen", 
                    "steelblue4", "lightblue", "goldenrod"),
                 c(summary(factor(taxanomy_filter_ordered$Order))))
  ) + scale_fill_manual(
    limits = as.factor(layout$name),
    values = rep(alpha(c( "black","yellow", "green", "gray", 
                          "darkgreen",  "pink", "hotpink",  "royalblue", 
                          "darkred", "firebrick2",  "bisque4", "orange" , 
                          "cyan3", "blue",  "darkolivegreen3", "plum1", 
                          "coral2", "purple",  "aquamarine3", "khaki4", 
                          "indianred", "mediumvioletred",  "tan4",
                          "turquoise4",  "wheat", "springgreen", 
                          "steelblue4", "lightblue", "goldenrod"),1),
                 c(summary(factor(taxanomy_filter_ordered$Order))))
  ) 
g1



########################
library(circlize)
names<-c( "o__","o__Acidaminococcales", "o__Bacteroidales",
          "o__Bifidobacteriales", "o__Burkholderiales", 
          "o__Christensenellales","o__Clostridia UCG-014", 
          "o__Clostridia vadinBB60 group"  , "o__Clostridiales",
          "o__Coriobacteriales" , "o__Corynebacteriales" ,
          "o__Desulfovibrionales" , "o__Enterobacterales" ,
          "o__Erysipelotrichales" , "o__Fusobacteriales" ,
          "o__Gastranaerophilales" , "o__Lachnospirales" ,
          "o__Lactobacillales", "o__Methanobacteriales" ,
          "o__Monoglobales" , "o__Oscillospirales" ,
          "o__Peptostreptococcales-Tissierellales" , "o__Pseudomonadales" ,
          "o__RF39" , "o__Rhizobiales",
          "o__Rhodospirillales" , "o__Veillonellales-Selenomonadales" ,
          "o__Verrucomicrobiales" , "o__Victivallales" )
colour<-c( "black","yellow", "green", "gray", 
           "darkgreen",  "pink", "hotpink",  "royalblue", 
           "darkred", "firebrick2",  "bisque4", "orange" , 
           "cyan3", "blue",  "darkolivegreen3", "plum1", 
           "coral2", "purple",  "aquamarine3", "khaki4", 
           "indianred", "mediumvioletred",  "tan4",
           "turquoise4",  "wheat", "springgreen", 
           "steelblue4", "lightblue", "goldenrod")
lgd = Legend(at = names, title = "", type = "points", 
             pch = 16, size = unit(5, "mm"), legend_gp = gpar(fill = colour, col = colour),ncol = 4,
             labels_gp = gpar(fontsize = 18), by_row=TRUE, background = "white")
draw(lgd)

#Legend
pdf("Figures/diet_legend_net.pdf",
    width=18, height=2)
op <- par(mar = rep(0, 4))
draw(lgd)
par(op)
dev.off()
