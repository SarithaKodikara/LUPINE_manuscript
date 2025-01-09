library(e1071)
library(ade4)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
set.seed(1234)
mantel.pvalue<-function(method="LUPINE",groups, minday=2, maxday=4){

  df= data.frame(days=rep(minday:maxday,2), rep(groups, each=maxday-1))
  mantel.pvalue_matrix <- matrix(0, dim(df)[1], dim(df)[1])

  # Calculate pairwise mantel.pvalue
  for (i in 1:(dim(df)[1]-1)) {
    for (j in (i+1):dim(df)[1]) {
      load(paste0("Results/",df[i,2],"/",method,"_Day",df[i,1],".rdata"))
      net_1<-(res$pvalue<0.05)*1
      net_1<-apply(net_1,c(1,2), function(x){ifelse(is.na(x),0,x)})
      lapl1<-as.dist(hamming.distance((net_1)))

      load(paste0("Results/",df[j,2],"/",method,"_Day",df[j,1],".rdata"))
      net_2<-(res$pvalue<0.05)*1
      net_2<-apply(net_2,c(1,2), function(x){ifelse(is.na(x),0,x)})
      lapl2<-as.dist(hamming.distance((net_2)))

      mantel.pvalue_matrix[i,j]<- mantel.pvalue_matrix[j,i]<-mantel.rtest(lapl1,lapl2)$pvalue
    }
  }
  return(mantel.pvalue_matrix)
}

group<-c("Plant", "Animal")
#res<-mantel.pvalue(method="LUPINE",groups= group, minday=2, maxday=15)
#save(res, file="Results/res_mantel.rdata")
load("Results/res_mantel.rdata")

rownames(res)<-rep(c(c("D[-3]","D[-2]","D[-1]"),paste0("D[",0:10,"]")),2)

combined_breaks <- c(seq(0, 0.05, 0.001), seq(0.0501,1,0.001))
combined_colors <- c(rev(colorRampPalette(brewer.pal(9,"RdPu")[1:5])(51)),
                     colorRampPalette(brewer.pal(9,"Blues"))(950))
combined_ramp <- colorRamp2(combined_breaks, combined_colors)
# Custom column annotation for label colors


col_ha <- HeatmapAnnotation(
  col_labels = anno_text(rownames(res),
                         rot = 90,
                         location=unit(1, 'npc'),
                         gp = gpar(col = c( rep('grey',3), rep('blue',5), rep('skyblue',6),
                                            rep('grey50',3) ,rep('tomato',5), rep('lightsalmon',6)),
                                   fontsize = 14)))

jpeg("Figures/FigC3.jpeg", units = 'in',
    width=10, height=10, res=300)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(res, col=combined_ramp,
                        show_heatmap_legend=F,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F,
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep('grey',3), rep('blue',5), rep('skyblue',6),
                                                     rep('grey50',3) ,rep('tomato',5), rep('lightsalmon',6))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", res[i, j]), x, y, gp = gpar(fontsize = 8))
                        })
par(op)
dev.off()
