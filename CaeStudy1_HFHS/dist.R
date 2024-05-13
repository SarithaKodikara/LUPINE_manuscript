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


method<-c( "LUPINE", "LUPINE_single")
group<-c("Normal","HFHS")

node_measures<-function(method="LUPINE",group,  day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("results/",group,"/",method,"_Day",day,".rdata"))
      net<-(res$pvalue<0.05)*1
      net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    net_ig<-graph_from_adjacency_matrix(net, mode = "undirected")
    return(net)
    
  }
  
}

g_ls.n<-lapply(2:4,function(k){node_measures("LUPINE","Normal",k)})
g_ls.h<-lapply(2:4,function(k){node_measures("LUPINE","HFHS",k)})
g_ls<-c(g_ls.n,g_ls.h)
dst<-nd.gdd(g_ls, out.dist = FALSE)$D
fit<-data.frame(cmdscale(dst, k=3))
names<-c(paste0("D[", c(1,4,7),"]"), paste0("D[", c(1,4,7),"]"))
fit$name<-c(names)
fit$color <- rep(c('#388ECC','#F68B33'), each=3)

g_ls.n1<-lapply(1:4,function(k){node_measures("LUPINE_single","Normal",k)})
g_ls.h1<-lapply(1:4,function(k){node_measures("LUPINE_single","HFHS",k)})
g_ls1<-c(g_ls.n1,g_ls.h1)
dst1<-nd.gdd(g_ls1, out.dist = FALSE)$D
fit1<-data.frame(cmdscale(dst1, k=3))
names<-c(paste0("D[", c(0,1,4,7),"]"), paste0("D[", c(0,1,4,7),"]"))
fit1$name<-c(names)
fit1$color <- rep(c('#388ECC','#F68B33'), each=4)


fit$title<-"LUPINE"
p1<-ggplot(fit, aes(x=X1, y=X2)) +
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+ 
  geom_point(size=3, col = fit$color)+ 
  labs(y = "MDS2", x = "MDS1")+xlim(c(-5.5,6.5))+ylim(-7,2.5)+
  theme_khakiwhite()+ facet_wrap(~title)

fit1$title<-"LUPINE_single"
p2<-ggplot(fit1, aes(x=X1, y=X2)) + 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "MDS2", x = "MDS1")+xlim(c(-4.5,7.5))+ylim(-3,7)+
  theme_khakiwhite()+ facet_wrap(~title)



pdf("Figures/HFHS_dis_LUPINE.pdf",
    width=3.5, height=3.5)
op <- par(mar = rep(0, 4))
p1
par(op)
dev.off()

pdf("Figures/HFHS_dis_LUPINE_single.pdf",
    width=3.5, height=3.5)
op <- par(mar = rep(0, 4))
p2
par(op)
dev.off()

#Legend
pdf("Figures/HFHS_legend.pdf",
    width=2.8, height=0.5)
op <- par(mar = rep(0, 4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Normal diet', 'HFHS diet'), pch=16, pt.cex=1.2, cex=1.2, 
       col = c('#388ECC','#F68B33'),box.lwd = 0.1,box.col = "grey",bg = "white", ncol=2 )
par(op)
dev.off()
