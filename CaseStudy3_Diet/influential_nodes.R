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

theme_bluewhite <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = alpha("lightblue",0.2)),
      panel.border = element_rect(color = "lightblue", fill = NA),
      axis.line = element_line(color = "lightblue"),
      axis.ticks = element_line(color = "lightblue"),
      axis.text = element_text(color = "black",family = "Courier",  size = (12)),
      plot.title = element_text(family = "Helvetica", face = "bold", size = (13)),
      axis.title = element_text(family = "Helvetica", size = (12)),
      strip.text = element_text(size=16),
      strip.text.x = element_text(margin = margin(0.3,0,0.3,0, "cm")),
      legend.key.size = unit(1.5, 'cm'),
      legend.key=element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=16)
    )
}

method<-c("LUPINE")
group<-c("Plant", "Animal")

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

number_of_edges<-function(method="LUPINE",group,  day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("Results/",group,"/",method,"_Day",day,".rdata"))
      net<-(res$pvalue<0.05)*1
      net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    return(sum(net)/2)
  }
}

df<-data.frame(Method=NA, Group=NA, Day=NA)
index=0
for(i in method){
  for(j in group){
    for(k in 1:15){
      if(!(k==1 &i=="LUPINE")){
        print(k)
        index=index+1
        df[index,1]<-i
        df[index,2]<-j
        df[index,3]<-k
        df[index,4:314]<-node_measures(i,j,k)
        
      }
    }
  }
}


#c( “#388ECC” “#F68B33" “#C2C2C2” “#009E73" “#CC79A7”)
bpls_res1<-pca(df[1:28,-(1:3)])
bpls_res1$names$sample<-c(c("D[-3]","D[-2]","D[-1]"),paste0("D[", 0:10,"]"), c("D[-3]","D[-2]","D[-1]"), paste0("D[", 0:10,"]"))
group<- c(rep(' Baseline',3), rep(' Plant',5), rep(' Plant_washout',6),
          rep(' Baseline',3), rep('Animal',5), rep('Animal_washout',6))
fig<-plotIndiv(bpls_res1,group=group, legend = F, 
          col.per.group = c('grey','blue',
                            'skyblue', 'tomato',
                            'lightsalmon'), cex = 5, title="LUPINE",
          xlim=c(-150,250))$graph

fit1<-data.frame(bpls_res1$variates$X)%>%cbind(name=bpls_res1$names$sample)
fit1$color<-c( rep('grey',3), rep('blue',5), rep('skyblue',6),
               rep('grey50',3) ,rep('tomato',5), rep('lightsalmon',6))
fit1$title = "LUPINE"
p1<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.2)+facet_wrap(~title)




pdf("Figures/Diet_ivi_LUPINE.1.pdf",
    width=5, height=5)
op <- par(mar = rep(0, 4))
p1
par(op)
dev.off()





load("data_modified/taxonomy_fil.rds")
load("data_modified/OTUdata_plant.rds")

Day0<-OTUdata_plant[,,1]
Colors=(brewer.pal(9,"RdPu"))
Colors=colorRampPalette(Colors)(99)

col1 = colorRamp2(seq(1,5,length=10), c("white", colorRampPalette(Colors)(9)))


#coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

col2 = list(Order = c( "o__" = "black","o__Acidaminococcales" = "yellow", "o__Bacteroidales" = "green",
                       "o__Bifidobacteriales" = "gray", "o__Burkholderiales" = "darkgreen", 
                       "o__Christensenellales" = "pink","o__Clostridia UCG-014" = "hotpink", 
                       "o__Clostridia vadinBB60 group" = "royalblue", "o__Clostridiales" = "darkred",
                       "o__Coriobacteriales" = "firebrick2", "o__Corynebacteriales" = "bisque4",
                       "o__Desulfovibrionales" = "orange" , "o__Enterobacterales" = "cyan3",
                       "o__Erysipelotrichales" = "blue", "o__Fusobacteriales" = "darkolivegreen3",
                       "o__Gastranaerophilales" = "plum1", "o__Lachnospirales" = "coral2",
                       "o__Lactobacillales" = "purple", "o__Methanobacteriales" = "aquamarine3",
                       "o__Monoglobales" = "khaki4", "o__Oscillospirales" = "indianred",
                       "o__Peptostreptococcales-Tissierellales" = "mediumvioletred", "o__Pseudomonadales" = "tan4",
                       "o__RF39" = "turquoise4", "o__Rhizobiales"= "wheat",
                       "o__Rhodospirillales" = "springgreen", "o__Veillonellales-Selenomonadales" = "steelblue4",
                       "o__Verrucomicrobiales" ="lightblue", "o__Victivallales" ="goldenrod"))
  


taxa_info<-taxanomy_filter_ordered[match(colnames(Day0), taxanomy_filter_ordered$OTU),]
taxa_info$Order<-factor(taxa_info$Order)
row_ha = rowAnnotation( Edges = anno_numeric(sapply(2:15, function(i){number_of_edges("LUPINE", "plant",i)}),
                                             bg_gp = gpar(fill = "orange", col = "red"),
                                         rg=c(0,2500)),
                        annotation_name_rot = 0
)
row_ha1 = rowAnnotation( Edges = anno_numeric(sapply(2:15, function(i){number_of_edges("LUPINE", "animal",i)}),
                                             bg_gp = gpar(fill = "orange", col = "red"),
                                             rg=c(0,2500)),
                        annotation_name_rot = 0
)

# Create the heatmap annotation
ha <- HeatmapAnnotation(
  Order = levels(taxa_info$Order),
  col = col2, show_legend = FALSE
)
m1<-as.matrix(df[1:14,-c(1:3)])
rownames(m1)<-c(paste0("P_D(", -3:-1, ")" ),paste0("P_D", 0:10))
colnames(m1)<-as.factor(rep(
  c( "black","yellow", "green", "gray", 
     "darkgreen",  "pink", "hotpink",  "royalblue", 
     "darkred", "firebrick2",  "bisque4", "orange" , 
     "cyan3", "blue",  "darkolivegreen3", "plum1", 
     "coral2", "purple",  "aquamarine3", "khaki4", 
     "indianred", "mediumvioletred",  "tan4",
     "turquoise4",  "wheat", "springgreen", 
     "steelblue4", "lightblue", "goldenrod"),
  c(summary(factor(taxanomy_filter_ordered$Order)))))

adjacency_df <- as.data.frame(as.table(m1))
# Use aggregate to sum values based on matching row and column names
sum_matrix <- aggregate(Freq ~ Var1 + Var2, data = adjacency_df, mean)

# Convert the result back to a matrix
sum_matrix <-as.matrix(log(reshape2::dcast(sum_matrix, Var1 ~ Var2, value.var = "Freq", fill = 0)[,-1]))
rownames(sum_matrix)<-c(paste0("P_D(", -3:-1, ")" ),paste0("P_D", 0:10))


m2<-as.matrix(df[15:28,-c(1:3)])
rownames(m2)<-c(paste0("A_D(", -3:-1, ")" ),paste0("A_D", 0:10))
colnames(m2)<-as.factor(rep(
  c( "black","yellow", "green", "gray", 
     "darkgreen",  "pink", "hotpink",  "royalblue", 
     "darkred", "firebrick2",  "bisque4", "orange" , 
     "cyan3", "blue",  "darkolivegreen3", "plum1", 
     "coral2", "purple",  "aquamarine3", "khaki4", 
     "indianred", "mediumvioletred",  "tan4",
     "turquoise4",  "wheat", "springgreen", 
     "steelblue4", "lightblue", "goldenrod"),
  c(summary(factor(taxanomy_filter_ordered$Order)))))

adjacency_df1 <- as.data.frame(as.table(m2))
# Use aggregate to sum values based on matching row and column names
sum_matrix1 <- aggregate(Freq ~ Var1 + Var2, data = adjacency_df1, mean)

# Convert the result back to a matrix
sum_matrix1 <-as.matrix(log(reshape2::dcast(sum_matrix1, Var1 ~ Var2, value.var = "Freq", fill = 0)[,-1]))
rownames(sum_matrix1)<-c(paste0("A_D(", -3:-1, ")" ),paste0("A_D", 0:10))

pdf("Figures/Diet_ivi_LUPINE.p.pdf",
    width=12, height=9)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(sum_matrix, cluster_rows = F, cluster_columns = F, col = col1,
        top_annotation = ha ,
        right_annotation = row_ha,
        column_split =levels(taxa_info$Order),
        row_split = c(rep('Baseline',3), rep('Diet',5), rep('Washout',6)),
        column_names_gp= gpar(fontsize = 0), 
        column_title_gp = gpar(fontsize = 10),
        show_heatmap_legend=FALSE,
        border=0, column_title_rot=90)

par(op)
dev.off()
pdf("Figures/Diet_ivi_LUPINE.a.pdf",
    width=12, height=9)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(sum_matrix1, cluster_rows = F, cluster_columns = F, col = col1,
                        top_annotation = ha ,
                        right_annotation = row_ha1,
                        column_split =levels(taxa_info$Order),
                        row_split =  c(rep('Baseline',3), rep('Diet',5), rep('Washout',6)),
                        column_names_gp= gpar(fontsize = 0), 
                        column_title_gp = gpar(fontsize = 10),
                        show_heatmap_legend=FALSE,
                        border=0, column_title_rot=90)

par(op)
dev.off()
