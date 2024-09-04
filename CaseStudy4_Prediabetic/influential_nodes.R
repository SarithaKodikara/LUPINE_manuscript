library(igraph)
library(mixOmics)
library(CINNA)
library(influential)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(ggplotify)
library(patchwork)
library(dplyr)

method<-c("LUPINE")
group<-c("PPT", "MDE")

node_measures<-function(method="LUPINE",group,  day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("Results/",group,"/",method,"_Day",day,".rdata"))
    net<-(res$pvalue<0.05)*1
    net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    net_ig<-graph_from_adjacency_matrix(net, mode = "undirected", weighted=NULL)
    IVI <- ivi(graph = net_ig)
    return(IVI)
  }
  
}


df<-data.frame(Method=NA, Group=NA, Day=NA)
index=0
for(i in method){
  for(j in group){
    k=2
      if(!(k==1 &i=="LUPINE")){
        print(k)
        index=index+1
        df[index,1]<-i
        df[index,2]<-j
        df[index,3]<-k
        df[index,4:94]<-node_measures(i,j,k)
    }
  }
}


load("data_modified/taxonomy.rds")
load("data_modified/data_PPT.rds")

Day0<-data_PPT[,,1]
Colors=(brewer.pal(9,"RdPu"))
Colors=colorRampPalette(Colors)(99)

col1 = colorRamp2(seq(1,5,length=10), c("white", colorRampPalette(Colors)(9)))


#coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

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
  


taxa_info<-taxonomy_order[match(colnames(Day0), taxonomy_order$Taxa),]
taxa_info$X6<-factor(taxa_info$X6)
# row_ha = rowAnnotation( Edges = anno_numeric(sapply(c("PPT","MDE"), function(i){number_of_edges("LUPINE", i,2)}),
#                                              bg_gp = gpar(fill = "orange", col = "red"),
#                                          rg=c(0,2500)),
#                         annotation_name_rot = 0
# )

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Genus = levels(taxa_info$X6),
  col = col2, show_legend = FALSE
)
m1<-as.matrix(df[1:2,-c(1:3)])
rownames(m1)<-c("PPT","MDE")
colnames(m1)<-as.factor(rep(
  c( "yellow", "green",  
     "darkgreen",  "pink", "hotpink",  "royalblue", 
     "darkred", "firebrick2",  "bisque4", "orange" , 
     "cyan3", "steelblue4",  "darkolivegreen3", "plum1", 
     "coral2", "purple",  "aquamarine3", "khaki4", 
     "indianred", "mediumvioletred",  "tan4",
     "turquoise4",  "wheat", "springgreen", "black"),
  c(summary(factor(taxa_info$X6)))))

adjacency_df <- data.frame(Taxa=rep(taxonomy_order$Taxa, each=2),as.table(m1))
adjacency_df <- adjacency_df %>% mutate(adjusted_value = ifelse(Var1 == "PPT", Freq, -Freq))

# Ensure the variable column is treated as a factor and retains its original order
adjacency_df$Taxa <- factor(adjacency_df$Taxa, levels = unique(adjacency_df$Taxa))

# Calculate y-position for annotations based on data
annotation_y <- length(unique(adjacency_df$Taxa)) + 4 # Adjust y position based on the number of variables


# Create the plot with colors from the 'color' column and unsorted y-axis
p<-ggplot(adjacency_df, aes(x = Taxa, y = adjusted_value, fill =  Var2)) +
  geom_bar(stat = "identity", position = "stack", width = 0.7)+
  coord_flip() +
  scale_y_continuous(breaks = seq(-100, 100, 10),
                     labels = abs(seq(-100, 100, 10))) +
  scale_x_discrete(limits = c(rep("",2),levels(adjacency_df$Taxa), rep("",2)))+# Show absolute values on the y-axis
  scale_fill_identity() +  # Use the actual colors from the 'color' column
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),  # Remove major horizontal grid lines
    panel.grid.minor.y = element_blank(),  # Remove minor horizontal grid lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.y=element_blank(),
    legend.position = "none"  # Remove the legend as colors are directly in the data
  ) +
  geom_hline(yintercept = 0, color="white", linewidth=2) +  # Add middle line at zero
  labs(x = "", y = "IVI")+
  # Add annotations for group labels
  annotate("text", x = annotation_y, y = 50, label = "PPT", 
           hjust = 0, fontface = "bold", color="grey33") +
  annotate("text", x = annotation_y, y = -50, label = "MDE", 
           hjust = 1, fontface = "bold", color="grey33")

pdf("Figures/Prediabetic_ivi_LUPINE.p.pdf",
    width=8, height=8)
op <- par(mar = rep(0, 4))
p
par(op)
dev.off()


