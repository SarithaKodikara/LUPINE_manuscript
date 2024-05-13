library(igraph)
library(mixOmics)
library(CINNA)
library(influential)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(tidyverse)
library(wesanderson)

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

method<-c("LUPINE", "LUPINE_single")
#method="sparcc"
group<-c("Normal", "HFHS")

node_measures<-function(method="LUPINE",group,  day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("Results/",group,"/",method,"_Day",day,".rdata"))
    net<-(res$pvalue<0.05)*1
    net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    net_ig<-graph_from_adjacency_matrix(net, mode = "undirected")
    IVI <- ivi(graph = net_ig)
    return(IVI)
  }
  
}

df<-data.frame(Method=NA, Group=NA, Day=NA)
index=0
for(i in method){
  for(j in group){
    for(k in 1:4){
      if(!(k==1 &i=="LUPINE")){
        print(k)
        index=index+1
        df[index,1]<-i
        df[index,2]<-j
        df[index,3]<-k
        df[index,4:215]<-node_measures(i,j,k)
        
      }
    }
  }
}

#c( “#388ECC” “#F68B33" “#C2C2C2” “#009E73" “#CC79A7”)
bpls_res1<-pca(df[1:6,-(1:3)])
bpls_res1$names$sample<-rep(c("D[1]", "D[4]", "D[7]"),2)
fig<-plotIndiv(bpls_res1,group=df$Group[1:6], legend = F, 
          col.per.group = c('#F68B33', '#388ECC'), cex = 5, title="LUPINE")$graph
fit1<-data.frame(bpls_res1$variates$X)%>%cbind(name=bpls_res1$names$sample)
fit1$color<-c(rep('#388ECC',3), rep('#F68B33',3))
fit1$title = "LUPINE"
p2<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-260,150)+xlim(-220, 260)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)



pca_res1<-pca(df[7:14,-(1:3)])
pca_res1$names$sample<-rep(c("D[0]", "D[1]", "D[4]", "D[7]"),2)
fig1<-plotIndiv(pca_res1,group=df$Group[7:14], legend = F, 
          col.per.group = c('#F68B33', '#388ECC'), cex = 5, title="LUPINE_single")$graph
fit1<-data.frame(pca_res1$variates$X)%>%cbind(name=pca_res1$names$sample)
fit1$color<-c(rep('#388ECC',4), rep('#F68B33',4))
fit1$title = "LUPINE_single"
p1<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-260,150)+xlim(-220, 260)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)



pdf("Figures/HFHS_ivi_LUPINE.pdf",
    width=3.5, height=3.5)
op <- par(mar = rep(0, 4))
p2
par(op)
dev.off()

pdf("Figures/HFHS_ivi_LUPINE_single.pdf",
    width=3.5, height=3.5)
op <- par(mar = rep(0, 4))
p1
par(op)
dev.off()


##### circle plot

data <- data.frame(
  taxa=paste( "Taxa ", seq(1,212), sep=""),
  group=factor(c(rep("green",54), rep("gray",2), rep("darkgreen",2),
                 rep("darkred",115), rep("firebrick2",7), rep("pink",1),
                 rep("tomato",1), rep("orange",3), rep("blue",15),
                 rep("purple",9), rep("hotpink",1), rep("lightblue",2)), 
               levels = c("green","gray","darkgreen",
                          "darkred","firebrick2","pink",
                          "tomato","orange","blue",
                          "purple","hotpink","lightblue"))
)
N_data<-cbind(data, t(df[1:3,4:215]))
colnames(N_data)[3:5]<-c("Day2", "Day4", "Day7")
H_data<-cbind(data, t(df[4:6,4:215]))
colnames(H_data)[3:5]<-c("Day2", "Day4", "Day7")

cPlotFunction<-function(data){
  # Transform data in a tidy format (long format)
  data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
  
  # Set a number of 'empty bar' to add at the end of each group
  empty_bar <- 2
  data$observation<-factor(data$observation,
                           levels = c("Day7", "Day4", "Day2"))
  nObsType <- nlevels(data$observation)
  to_add <- data.frame( matrix(NA, empty_bar*nlevels(data$group)*nObsType, ncol(data)) )
  colnames(to_add) <- colnames(data)
  to_add$group <- rep(levels(data$group), each=empty_bar*nObsType )
  data <- rbind(data, to_add)
  data <- data %>% arrange(group, taxa)
  data$id <- rep( seq(1, nrow(data)/nObsType) , each=nObsType)
  
  # Get the name and the y position of each label
  label_data <- data %>% group_by(id, taxa) %>% summarize(tot=sum(value))
  number_of_bar <- nrow(label_data)
  angle <- 90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)
  label_data$hjust <- ifelse( angle < -90, 1, 0)
  label_data$angle <- ifelse(angle < -90, angle+180, angle)
  
  # prepare a data frame for base lines
  base_data <- data %>% 
    group_by(group) %>% 
    summarize(start=min(id), end=max(id) - empty_bar) %>% 
    rowwise() %>% 
    mutate(title=mean(c(start, end)))
  
  # prepare a data frame for grid (scales)
  grid_data<-data.frame(start=rep(NA,59),end=rep(NA,59))
  grid_data$end <- seq(1,234,4)-1
  grid_data$start <- seq(1,234,4)
  
  # Make the plot
  p <- ggplot(data) +      
    
    # Add the stacked bar
    geom_bar(aes(x=as.factor(id), y=value, fill=observation), stat="identity", alpha=0.8) +
    scale_fill_manual(values = wes_palette("GrandBudapest2", n = 3))+
    
    # Add a val=100/75/50/25 lines. I do it at the beginning to make sur barplots are OVER it.
    #geom_segment(data=grid_data, aes(x = end, y = 0, xend = start, yend = 0), colour = "grey", alpha=1, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 50, xend = start, yend = 50), colour = "black", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = end, y = 100, xend = start, yend = 100), colour = "black", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = start, y = 150, xend = end, yend = 150), colour = "black", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    geom_segment(data=grid_data, aes(x = start, y = 200, xend = end, yend = 200), colour = "black", alpha=0.5, size=0.3 , inherit.aes = FALSE ) +
    
    # Add text showing the value of each 100/75/50/25 lines
    ggplot2::annotate("text", x = rep(max(data$id),6), y = c(0, 50, 100, 150, 200, 250), label = c("0", "50", "100", "150", "200", 250) , color="grey", size=5 , angle=0, fontface="bold", hjust=0.9) +
    
    ylim(-200,200) +
    theme_minimal() +
    theme(
      legend.position = "none",
      axis.text = element_blank(),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      plot.margin = unit(rep(-1,4), "cm") 
    ) +
    coord_polar() +
    
    # Add labels on top of each bar
    #geom_text(data=label_data, aes(x=id, y=tot+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=5, angle= label_data$angle, inherit.aes = FALSE ) +
    
    # Add base line information
    geom_segment(data=base_data, aes(x = start, y = -5, xend = end, yend = -5), 
                 colour = c("green","gray","darkgreen",
                            "darkred","firebrick2","pink",
                            "tomato","orange","blue",
                            "purple","hotpink","lightblue"), 
                 alpha=1, size=2 , inherit.aes = FALSE )  
  
  
  return(p)
}


pdf("Figures/HFHS_cir_LUPINE_normal.pdf",
    width=5, height=5)
op <- par(mar = rep(0, 4))
cPlotFunction(N_data)
par(op)
dev.off()


pdf("Figures/HFHS_cir_LUPINE_hfhs.pdf",
    width=5, height=5)
op <- par(mar = rep(0, 4))
cPlotFunction(H_data)
par(op)
dev.off()

#Legend
pdf("Figures/BEME_legend_cir.pdf",
    width=0.9, height=0.9)
op <- par(mar = rep(0, 4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Day 1', 'Day 4', 'Day 7'), pch=15, pt.cex=1.2, cex=1.2, 
       col = c( "#D8A499","#C6CDF7", "#E6A0C4"),box.lwd = 0.1,box.col = "grey",bg = "white" )
par(op)
dev.off()
