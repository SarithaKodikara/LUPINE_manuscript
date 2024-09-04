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

set.seed(123456)
load("data_modified/taxonomy.rds")

load("Results/MDE/LUPINE_Day2.rdata")
netMDE<-(res$pvalue<0.05)*1
netMDE<-apply(netMDE,c(1,2), function(x){ifelse(is.na(x),0,x)})

load("Results/PPT/LUPINE_Day2.rdata")
netPPT<-(res$pvalue<0.05)*1
netPPT<-apply(netPPT,c(1,2), function(x){ifelse(is.na(x),0,x)})

net=netPPT

g <- graph_from_adjacency_matrix(net, mode="undirected", weighted=NULL)

# provide some names
V(g)$name <- 1:vcount(g)



# plot using ggraph
graph_tbl <- g %>%
  as_tbl_graph() %>%
  activate(nodes) %>%
  mutate(degree  = centrality_degree()) %>%
  mutate(community = as.factor(rep(
    c( "yellow", "green",  
       "darkgreen",  "pink", "hotpink",  "royalblue", 
       "darkred", "firebrick2",  "bisque4", "orange" , 
       "cyan3", "steelblue4",  "darkolivegreen3", "plum1", 
       "coral2", "purple",  "aquamarine3", "khaki4", 
       "indianred", "mediumvioletred",  "tan4",
       "turquoise4",  "wheat", "springgreen", "black"),
    c(summary(factor(taxonomy_order$X6))))))


##############
layout1<-as.matrix(create_layout(induced_subgraph(graph_tbl, 1:3), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout2<-as.matrix(create_layout(induced_subgraph(graph_tbl, 4), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2) 
layout3<-as.matrix(create_layout(induced_subgraph(graph_tbl,5:14), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout4<-as.matrix(create_layout(induced_subgraph(graph_tbl,15:16), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout5<-as.matrix(create_layout(induced_subgraph(graph_tbl,17), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout6<-as.matrix(create_layout(induced_subgraph(graph_tbl,18:28), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout7<-as.matrix(create_layout(induced_subgraph(graph_tbl,29), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout8<-as.matrix(create_layout(induced_subgraph(graph_tbl,30), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout9<-as.matrix(create_layout(induced_subgraph(graph_tbl,31), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout10<-as.matrix(create_layout(induced_subgraph(graph_tbl,32), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout11<-as.matrix(create_layout(induced_subgraph(graph_tbl,33:42), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout12<-as.matrix(create_layout(induced_subgraph(graph_tbl,43:44), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout13<-as.matrix(create_layout(induced_subgraph(graph_tbl,45:46), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout14<-as.matrix(create_layout(induced_subgraph(graph_tbl,47:51), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout15<-as.matrix(create_layout(induced_subgraph(graph_tbl,52:62), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout16<-as.matrix(create_layout(induced_subgraph(graph_tbl,63:65), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout17<-as.matrix(create_layout(induced_subgraph(graph_tbl,66), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout18<-as.matrix(create_layout(induced_subgraph(graph_tbl,67), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout19<-as.matrix(create_layout(induced_subgraph(graph_tbl,68), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout20<-as.matrix(create_layout(induced_subgraph(graph_tbl,69), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout21<-as.matrix(create_layout(induced_subgraph(graph_tbl,70), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout22<-as.matrix(create_layout(induced_subgraph(graph_tbl,71:72), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout23<-as.matrix(create_layout(induced_subgraph(graph_tbl,73:76), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout24<-as.matrix(create_layout(induced_subgraph(graph_tbl,77), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)
layout25<-as.matrix(create_layout(induced_subgraph(graph_tbl,78:91), layout = 'igraph', algorithm = 'sphere')[,1:2], ncol=2)

layout0 <- list(layout1, layout2, layout3, layout4,
              layout5,layout6,layout7,layout8,
              layout9,layout10,layout11,layout12,
              layout13,layout14,layout15,layout16,
              layout17,layout18,layout19,layout20,
              layout21,layout22,layout23,layout24,
              layout25)

graphs<-list(induced_subgraph(graph_tbl, 1:3),
             induced_subgraph(graph_tbl, 4),
             induced_subgraph(graph_tbl, 5:14),
             induced_subgraph(graph_tbl, 15:16),
             induced_subgraph(graph_tbl, 17),
             induced_subgraph(graph_tbl, 18:28),
             induced_subgraph(graph_tbl, 29),
             induced_subgraph(graph_tbl,30),induced_subgraph(graph_tbl,31),
             induced_subgraph(graph_tbl,32),induced_subgraph(graph_tbl,33:42),
             induced_subgraph(graph_tbl,43:44),induced_subgraph(graph_tbl,45:46),
             induced_subgraph(graph_tbl,47:51),induced_subgraph(graph_tbl,52:62),
             induced_subgraph(graph_tbl,63:65),induced_subgraph(graph_tbl,66),
             induced_subgraph(graph_tbl,67),induced_subgraph(graph_tbl,68),
             induced_subgraph(graph_tbl,69),induced_subgraph(graph_tbl,70),
             induced_subgraph(graph_tbl,71:72),induced_subgraph(graph_tbl,73:76),
             induced_subgraph(graph_tbl,77),induced_subgraph(graph_tbl,78:91))
lay <- merge_coords(graphs, layout0)

layout <- create_layout(graph_tbl, layout = 'igraph', algorithm = 'circle')
layout$x<-lay[,1]+layout$x*30
layout$y<-lay[,2]+layout$y*20

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
    values = rep(c( "yellow", "green",  
                    "darkgreen",  "pink", "hotpink",  "royalblue", 
                    "darkred", "firebrick2",  "bisque4", "orange" , 
                    "cyan3", "steelblue4",  "darkolivegreen3", "plum1", 
                    "coral2", "purple",  "aquamarine3", "khaki4", 
                    "indianred", "mediumvioletred",  "tan4",
                    "turquoise4",  "wheat", "springgreen","black"),
                 c(summary(factor(taxonomy_order$X6))))
  ) + scale_fill_manual(
    limits = as.factor(layout$name),
    values = rep(alpha(c( "yellow", "green",
                          "darkgreen",  "pink", "hotpink",  "royalblue", 
                          "darkred", "firebrick2",  "bisque4", "orange" , 
                          "cyan3", "steelblue4",  "darkolivegreen3", "plum1", 
                          "coral2", "purple",  "aquamarine3", "khaki4", 
                          "indianred", "mediumvioletred",  "tan4",
                          "turquoise4",  "wheat", "springgreen","black"),1),
                 c(summary(factor(taxonomy_order$X6))))
  ) 
g1

pdf("Figures/PPT.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
g1

par(op)
dev.off()



########################
library(circlize)
names<-levels(factor(taxonomy_order$X6))
colour<-c( "yellow", "green",
           "darkgreen",  "pink", "hotpink",  "royalblue", 
           "darkred", "firebrick2",  "bisque4", "orange" , 
           "cyan3", "steelblue4",  "darkolivegreen3", "plum1", 
           "coral2", "purple",  "aquamarine3", "khaki4", 
           "indianred", "mediumvioletred",  "tan4",
           "turquoise4",  "wheat", "springgreen","black")
lgd = Legend(at = names, title = "", type = "points", 
             pch = 16, size = unit(5, "mm"), legend_gp = gpar(fill = colour, col = colour),ncol = 4,
             labels_gp = gpar(fontsize = 18), by_row=TRUE, background = "white")
#draw(lgd)

#Legend
pdf("Figures/prediabetic_legend_net.pdf",
    width=14, height=2)
op <- par(mar = rep(0, 4))
draw(lgd)
par(op)
dev.off()


## Aggregate 

df <- as.data.frame(as.table(netMDE-netPPT))
# Rename columns for clarity
colnames(df) <- c('From', 'To', 'Weight')
geneus <- c(taxonomy_order$X6)

# Add factors
df$From_Group <- geneus[df$From]
df$To_Group <- geneus[df$To]

# Aggregate weights based on group pairs
agg_df <- df %>%
  group_by(From_Group, To_Group) %>%
  summarise(Aggregated_Weight = sum(Weight)) %>%
  ungroup()
# Create a matrix of unique groups
groups <- unique(c(agg_df$From_Group, agg_df$To_Group))
agg_matrix <- matrix(0, nrow = length(groups), ncol = length(groups))
rownames(agg_matrix) <- colnames(agg_matrix) <- groups

# Fill the aggregated matrix
for (i in seq_len(nrow(agg_df))) {
  from <- agg_df$From_Group[i]
  to <- agg_df$To_Group[i]
  agg_matrix[from, to] <- agg_matrix[from, to] + agg_df$Aggregated_Weight[i]
  agg_matrix[to, from] <- agg_matrix[to, from] + agg_df$Aggregated_Weight[i] # Symmetric addition
}

print(agg_matrix/2)

Colors=(brewer.pal(11,"PRGn"))
Colors=colorRampPalette(Colors)(100)

col1 = colorRamp2(seq(-20,20,length=100), Colors)
col2 = list(Genus = c( "g__Alistipes" = "yellow", "g__Anaerostipes" = "green", "g__Bacteroides" = "darkgreen", 
                       "g__Bifidobacterium" = "pink","g__Bilophila" = "hotpink", 
                       "g__Blautia" = "royalblue", "g__Burkholderiales_unclassified" = "darkred",
                       "g__Butyricicoccus" = "firebrick2", "g__Butyricimonas" = "bisque4",
                       "g__Clostridiales_unclassified" = "orange" , "g__Clostridium" = "cyan3",
                       "g__Coprococcus" = "steelblue4", "g__Dorea" = "darkolivegreen3",
                       "g__Eubacterium" = "plum1", "g__Faecalibacterium" = "coral2",
                       "g__Firmicutes_unclassified" = "purple", "g__Fusicatenibacter" = "aquamarine3",
                       "g__Gemmiger" = "khaki4", "g__Lachnospiraceae_unclassified" = "indianred",
                       "g__Odoribacter" = "mediumvioletred", "g__Oscillibacter" = "tan4",
                       "g__Parabacteroides" = "turquoise4", "g__Roseburia"= "wheat",
                       "g__Ruminococcus" = "springgreen", "g__unknown"="black"))


ha <- HeatmapAnnotation(
  Genus = levels(factor(taxonomy_order$X6)),
  col = col2 , show_legend = FALSE,
  annotation_name_gp =gpar(fontsize = 0)
)
row_ha <- rowAnnotation(
  Genus = levels(factor(taxonomy_order$X6)),
  col = col2 , show_legend = FALSE,
  annotation_name_gp =gpar(fontsize = 0)
)



pdf("Figures/Prediabetic_connection_diff.pdf",
    width=12, height=9)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(agg_matrix/2, col=col1,
                        column_names_gp= gpar(fontsize = 10), 
                        column_title_gp = gpar(fontsize = 10),
                        show_heatmap_legend=T,
                        border=0, 
                        top_annotation = ha,
                        right_annotation = row_ha, 
                        cluster_rows = T, cluster_columns = T)

par(op)
dev.off()
