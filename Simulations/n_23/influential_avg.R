library(igraph)
library(mixOmics)
library(CINNA)
library(influential)
library(ggplot2)
library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
library(patchwork)
library(ggplotify)
library(abind)

method<-c( "bPLS", "PCA", "SpiecEasi_mb", "SpiecEasi_glasso", "sparcc")
para<-readRDS("Sim_para/para.rds")

#True1<-ivi(graph = graph_from_adjacency_matrix(para$binary.graph1_n, 
                                               #mode = "undirected")) 
#True2<-ivi(graph = graph_from_adjacency_matrix(para$binary.graph1_h, 
                                              # mode = "undirected"))

node_measures<-function(method="bPLS",iter,  day){
  if(!(day==1 &method=="bPLS")){
    load(paste0("Results/",method,"/iter",iter,"_d",day,".rdata"))

    if(method%in%c("SpiecEasi_mb","SpiecEasi_glasso")){
      net<-res$Graph
    }else{
      #cutoff.value<-cutoff(day,res$pvalue)
      net<-(res$pvalue>0.95)*1
      net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})

    }
    net_ig<-graph_from_adjacency_matrix(net,
                                        mode = "undirected")
    IVI <-ivi(graph = net_ig)
    return(IVI)
  }

}

influentialNodes<-function(method){
  g_ls<-lapply(1:50, function(iter){
    sapply(2:10,function(k){print(paste(iter,"_",k));node_measures(method,iter,k)})}
    )
  abind_ls <- do.call(abind, c(g_ls, list(along=3)))
  df<-t(apply(abind_ls, 1:2, mean))
  return(df)
}

group<-c(rep(c('B','A'), each=4),"A","T1","T2")

# fit_bPLS<-influentialNodes("bPLS")
# fit_PCA<-influentialNodes("PCA")
# fit_sparcc<-influentialNodes("sparcc")
# fit_mb<-influentialNodes("SpiecEasi_mb")
# fit_glasso<-influentialNodes("SpiecEasi_glasso")
#res_lst<-list(fit_bPLS, fit_PCA,fit_sparcc,fit_mb,fit_glasso, True1, True2)
#save(res_lst, file="Results/Results_inf_avg.rdata")
load("Results/Results_inf_avg.rdata")
fit_bPLS<-res_lst[[1]]
fit_PCA<-res_lst[[2]]
fit_sparcc<-res_lst[[3]]
fit_mb<-res_lst[[4]]
fit_glasso<-res_lst[[5]]
True1<-res_lst[[6]]
True2<-res_lst[[7]]



mynamestheme <- theme(
  plot.title = element_text(family = "Helvetica", face = "bold", size = (18)),
  axis.title = element_text(family = "Helvetica", size = (12)),
  axis.text = element_text(family = "Courier",  size = (12))
)

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

pca_res1<-pca(rbind(fit_PCA, True1, True2))
pca_res1$names$sample<-c(paste0("D[", 2:10,"]"),"T[1]","T[2]")
fig<-plotIndiv(pca_res1, group=group, legend = F,  
               col.per.group =c('#F68B33', '#388ECC','#C2C2C2', '#009E73'), 
               cex = 5, title="LUPINE_single",
               ylim=c(-105,100), xlim = c(-155, 155))$graph

fit1<-data.frame(pca_res1$variates$X)%>%cbind(name=pca_res1$names$sample)
fit1$color<-c(rep('#F68B33',4), rep('#388ECC',5),'#C2C2C2', '#009E73')
fit1$title = "LUPINE_single"
p1<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-140,100)+xlim(-155, 165)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)



bpls_res1<-pca(rbind(fit_bPLS, True1, True2))
bpls_res1$names$sample<-c(paste0("D[", 2:10,"]"),"T[1]","T[2]")
fig<-plotIndiv(bpls_res1, group=group, legend = F,  
                         col.per.group =c('#F68B33', '#388ECC','#C2C2C2', '#009E73'), 
                         cex = 5, title="LUPINE",ylim=c(-105,100), 
               xlim = c(-155, 165))$graph
fit1<-data.frame(bpls_res1$variates$X)%>%cbind(name=bpls_res1$names$sample)
fit1$color<-c(rep('#F68B33',4), rep('#388ECC',5),'#C2C2C2', '#009E73')
fit1$title = "LUPINE"
p2<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-140,100)+xlim(-155, 155)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)



sparcc_res1<-pca(rbind(fit_sparcc, True1, True2))
sparcc_res1$names$sample<-c(paste0("D[", 2:10,"]"),"T[1]","T[2]")
fig<-plotIndiv(sparcc_res1, group=group, legend = F,  
               col.per.group =c('#F68B33', '#388ECC','#C2C2C2', '#009E73'), 
               cex = 5, title="sparcc",ylim=c(-105,100),
               xlim = c(-155, 165))$graph
fit1<-data.frame(sparcc_res1$variates$X)%>%cbind(name=sparcc_res1$names$sample)
fit1$color<-c(rep('#F68B33',4), rep('#388ECC',5),'#C2C2C2', '#009E73')
fit1$title = "sparcc"
p3<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-140,100)+xlim(-155, 155)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)


mb_res1<-pca(rbind(fit_mb, True1, True2))
mb_res1$names$sample<-c(paste0("D[", 2:10,"]"),"T[1]","T[2]")
fig<-plotIndiv(mb_res1, group=group, legend = F,  
               col.per.group =c('#F68B33', '#388ECC','#C2C2C2', '#009E73'), 
               cex = 5, title="SpiecEasi_mb",ylim=c(-105,100), 
               xlim = c(-155, 165))$graph
fit1<-data.frame(mb_res1$variates$X)%>%cbind(name=mb_res1$names$sample)
fit1$color<-c(rep('#F68B33',4), rep('#388ECC',5),'#C2C2C2', '#009E73')
fit1$title = "SpiecEasi_mb"
p4<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-140,100)+xlim(-155, 155)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)

glasso_res1<-pca(rbind(fit_glasso, True1, True2))
glasso_res1$names$sample<-c(paste0("D[", 2:10,"]"),"T[1]","T[2]")
fig<-plotIndiv(glasso_res1, group=group, legend = F,  
               col.per.group =c('#F68B33', '#388ECC','#C2C2C2', '#009E73'), 
               cex = 5, title="SpiecEasi_glasso",ylim=c(-105,100), xlim = c(-155, 155))$graph
fit1<-data.frame(glasso_res1$variates$X)%>%cbind(name=glasso_res1$names$sample)
fit1$color<-c(rep('#F68B33',4), rep('#388ECC',5),'#C2C2C2', '#009E73')
fit1$title = "SpiecEasi_glasso"
p5<-ggplot(fit1, aes(x=PC1, y=PC2)) + 
  ylim(-140,100)+xlim(-155, 165)+
  geom_point(size=3, col = fit1$color)+ 
  labs(y = "PC1", x = "PC2")+
  theme_bluewhite()+ 
  geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5)+facet_wrap(~title)


pdf("Figures/n23_pca1_n.pdf",
    width=3, height=3)
set.seed(12345)
p1
dev.off()

Colors=(brewer.pal(9,"RdPu"))
Colors=colorRampPalette(Colors)(100)
col = colorRamp2(seq(1,100,length=100), colorRampPalette(Colors)(100))



m1<-(rbind(True1,fit_sparcc[1:4,],fit_mb[1:4,], fit_glasso[1:4,],
           fit_PCA[1:4,], fit_bPLS[1:4,]))
rownames(m1)<-c(" ",paste0("D",2:5),paste0("D",2:5),
                paste0("D",2:5),paste0("D",2:5),paste0("D",2:5))

pdf("Figures/sim_n23_netdis2.2_n.pdf",
    width=8, height=5)
op <- par(mar = rep(0, 4))

rowSplit<-factor(c("True",rep("SparCC",4),rep("SpiecEasi_mb",4), 
               rep("SpiecEasi_glasso",4), rep("LUPINE_single",4), rep("LUPINE",4)),
             levels = c("True","SparCC","SpiecEasi_mb", "SpiecEasi_glasso","LUPINE_single", "LUPINE"))

Heatmap(m1, cluster_rows = F, cluster_columns = F, col = col,
        row_split = rowSplit,
        column_names_gp= gpar(fontsize = 8), 
        row_title_gp = gpar(fontsize = 18),
        #row_title =NULL, 
        row_title_rot = 0,
        show_heatmap_legend=FALSE,
        border=1, column_names_rot=90,
        column_title = "Nodes influence for days 2 to 5")
par(op)
dev.off()



m1<-(rbind(True2,fit_sparcc[5:9,],fit_mb[5:9,], fit_glasso[5:9,],
           fit_PCA[5:9,], fit_bPLS[5:9,]))
rownames(m1)<-c(" ",paste0("D",6:10),paste0("D",6:10),
                paste0("D",6:10),paste0("D",6:10),paste0("D",6:10))

pdf("Figures/sim_n23_netdis2.3_n.pdf",
    width=8, height=5)
op <- par(mar = rep(0, 4))
rowSplit<-factor(c("True",rep("SparCC",5),rep("SpiecEasi_mb",5), 
                   rep("SpiecEasi_glasso",5), rep("LUPINE_single",5), rep("LUPINE",5)),
                 levels = c("True","SparCC","SpiecEasi_mb", "SpiecEasi_glasso","LUPINE_single", "LUPINE"))


Heatmap(m1, cluster_rows = F, cluster_columns = F, col = col,
        row_split = rowSplit,
        column_names_gp= gpar(fontsize = 8), 
        row_title_gp = gpar(fontsize = 18),
        #row_title =NULL, 
        row_title_rot = 0,
        show_heatmap_legend=FALSE,
        border=1, column_names_rot=90,
        column_title = "Nodes influence for days 6 to 10")
par(op)
dev.off()

