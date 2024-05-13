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
library(tidyverse)
library(wesanderson)

method<-c("LUPINE")

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

node_measures<-function(method="LUPINE", day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("Results/",method,"_Day",day,".rdata"))
    net<-(res$pvalue<0.05)*1
    net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    net_ig<-graph_from_adjacency_matrix(net, mode = "undirected")
    IVI <- ivi(graph = net_ig)
    return(IVI)
  }
  
}

df<-data.frame(Method=NA, Day=NA)
index=0
for(i in method){
    for(k in 1:10){
      if(!(k==1 &i=="LUPINE")){
        print(k)
        index=index+1
        df[index,1]<-i
        df[index,2]<-k
        df[index,3:128]<-node_measures(i,k)
        
      }
    }
}

#c( “#388ECC” “#F68B33" “#C2C2C2” “#009E73" “#CC79A7”)
bpls_res1<-pca(df[1:9,-(1:2)])
bpls_res1$names$sample<-paste0("D[", c(1:2,5:7, 9, 12:14), "]")
group<- c(rep(' naive', 3),rep('anti', 2),rep('VRE', 4))
fig<-plotIndiv(bpls_res1,group=group, legend = F, 
          col.per.group = c('#009E73','#F68B33', '#388ECC'), 
          cex = 5,title="LUPINE")$graph
fit1<-data.frame(bpls_res1$variates$X)%>%cbind(name=bpls_res1$names$sample)
fit1$color<-c(rep('#009E73',3), rep('#F68B33',2), rep('#388ECC',4))
fit1$title = "LUPINE"
p1<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-120,150)+xlim(-170, 120)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)


pdf("Figures/VRE_ivi_LUPINE.1.pdf",
    width=5, height=5)
op <- par(mar = rep(0, 4))
p1
par(op)
dev.off()


##### Line plots

load("data_filtered/taxonomy.rds")
data <- data.frame(
  taxa=paste( "Taxa ", seq(1,126), sep=""),
  group=factor(rep(levels(factor(taxanomy_filter_ordered$V6)),
                   c(summary(factor(taxanomy_filter_ordered$V6)))))
)
data<-cbind(data, t(df[1:9,3:128]))
colnames(data)[3:11]<-c("D_1", "D_2", "D_5", "D_6",
                        "D_7", "D_9", "D_12", "D_13","D_14")

# Transform data in a tidy format (long format)
data <- data %>% gather(key = "observation", value="value", -c(1,2)) 
data$day<-as.integer(str_replace(data$observation,"D_", ""))

ggplot(data,
       aes(x = day, y = value,
           group = taxa, col=group))+
  geom_line()+ facet_wrap(. ~ group, ncol=3)+
  scale_colour_manual(values = c("black","grey48","pink2", "green", "darkred", "orange","red","deepskyblue",
                                 "purple", "hotpink","darkgreen","yellow","tomato", "blue")) +
  theme(legend.position = 'none')


pdf("Figures/ivi_LUPINE_line.pdf",
    width=5, height=5)
ggplot(data) +
  geom_rect(data=data.frame(variable=factor(1:3)),
            aes(xmin = 0,xmax = 5,ymin = -Inf, ymax = Inf),
            fill='#009E73', alpha = .1)+ylab('IVI')+
  geom_rect(data=data.frame(variable=factor(1:3)),
            aes(xmin = 5,xmax = 7,ymin = -Inf, ymax = Inf),
            fill='#F68B33', alpha = .1) +
  geom_rect(data=data.frame(variable=factor(1:3)),
            aes(xmin = 7,xmax = 14,ymin = -Inf, ymax = Inf),
            fill='#388ECC', alpha = .1) +
  geom_line(aes(x = day, y = value,
                group = taxa, col=group)) + facet_wrap(. ~ group, ncol=3)+
  scale_colour_manual(values = c("black","grey48","pink2", "green", "darkred", "orange","red","deepskyblue",
                                 "purple", "hotpink","darkgreen","yellow","tomato", "blue")) +
  theme(legend.position = 'none')
dev.off()





