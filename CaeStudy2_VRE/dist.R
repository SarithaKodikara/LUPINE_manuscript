library(igraph)
library(mixOmics)
library(readr)
library(NetworkDistance)
library(patchwork)

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

method<-c( "LUPINE")

node_measures<-function(method="LUPINE",day){
  if(!(day==1 &method=="LUPINE")){
    load(paste0("Results/",method,"_Day",day,".rdata"))
    net<-(res$pvalue<0.05)*1
    net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
    return(net)
  }
}

g_ls<-lapply(2:10,function(k){node_measures("LUPINE",k)})
dst<-nd.gdd(g_ls, out.dist = FALSE)$D
fit<-data.frame(cmdscale(dst, k=3))
names<-paste0("D[", c(1:2,5:7, 9, 12:14),"]")
fit$name<-factor(names, level=names)
fit$phase<-c(rep(' Naive', 3),rep('Antibiotic', 2),rep('VRE', 4))


fit$title = "LUPINE"
p1<-ggplot(fit, aes(x=X1, y=X2, color=phase)) + 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+ 
  geom_point(size=3)+ 
  ylim(-2.5,3.5)+xlim(-5.5, 3)+ 
  labs(y = "MDS2", x = "MDS1")+
  theme_khakiwhite()+ 
  scale_colour_manual(values = c('#009E73', '#F68B33', '#388ECC'))+facet_wrap(~title)+
  theme(legend.position ='none' )
  


pdf("Figures/VRE_dis_LUPINE.pdf",
    width=5, height=5)
op <- par(mar = rep(0, 4))
p1
par(op)
dev.off()


#Legend
pdf("Figures/VRE_legend.pdf",
    width=1.2, height=1)
op <- par(mar = rep(0, 4))
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("topleft", legend =c('Naive', 'Antibiotic', 'VRE'), pch=16, pt.cex=1.2, cex=1.2, 
       col = c('#009E73', '#F68B33', '#388ECC'),box.lwd = 0.1,box.col = "grey",bg = "white", ncol=1 )
par(op)
dev.off()
