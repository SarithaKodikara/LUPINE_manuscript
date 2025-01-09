library(igraph)
library(mixOmics)
library(readr)
library(NetworkDistance)


theme_khakiwhite <- function (base_size = 11, base_family = "") {
  theme_bw() %+replace% 
    theme(
      panel.grid.major  = element_line(color = "white"),
      panel.background = element_rect(fill = alpha("khaki2",0.2)),
      panel.border = element_rect(color = "khaki2", fill = NA),
      axis.line = element_line(color = "khaki2"),
      axis.ticks = element_line(color = "khaki2"),
      axis.text = element_text(color = "black",family = "Courier",  size = (12)),
      plot.title = element_text(family = "Helvetica", face = "bold", size = (18)),
      axis.title = element_text(family = "Helvetica", size = (12)),
      strip.text = element_text(size=16),
      strip.text.x = element_text(margin = margin(0.3,0,0.3,0, "cm")),
      legend.key.size = unit(1.5, 'cm'),
      legend.key=element_blank(),
      legend.title = element_blank(),
      legend.text = element_text(size=16)
    )
}


node_measures<-function(method="LUPINE",group,  day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("Results/",group,"/",method,"_Day",day,".rdata"))
    net<-(res$pvalue<0.05)*1
    net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    net_ig<-graph_from_adjacency_matrix(net, mode = "undirected")
    return(net)
    
  }
  
}

g_ls.p<-lapply(1:2,function(k){node_measures("LUPINEsin","PPT",k)})
g_ls.m<-lapply(1:2,function(k){node_measures("LUPINEsin","MDE",k)})
g_ls<-c(g_ls.p,g_ls.m)
dst<-nd.gdd(g_ls, out.dist = FALSE)$D


fit<-data.frame(cmdscale(dst, k=3))
names<-rep(paste0("D[", 0:1,"]"),2)
fit$name<-c(names)
fit$color <- c(  'skyblue','blue', 'lightsalmon','tomato')


fit$title<-"LUPINE"
p1<-ggplot(fit, aes(x=X1, y=X2)) + 
  ylim(-3,4)+
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
  geom_point(size=3, col = fit$color)+ 
  labs(y = "MDS2", x = "MDS1")+
  theme_khakiwhite()+ facet_wrap(~title)+ 
  labs(colour = "Diet_phase")+
  theme(legend.position ='none' )
p1



