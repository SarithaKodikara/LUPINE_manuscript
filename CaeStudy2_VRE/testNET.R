library(e1071)
library(ade4)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
set.seed(1234)
mantel.pvalue<-function(method="LUPINE",groups, minday=2, maxday=4){

  df= data.frame(days=c(2:maxday))
  mantel.pvalue_matrix <- matrix(0, dim(df)[1], dim(df)[1])

  # Calculate pairwise mantel.pvalue
  for (i in 1:(dim(df)[1]-1)) {
    for (j in (i+1):dim(df)[1]) {

      load(paste0("Results/",method,"_Day",df[i,1],".rdata"))
      net_1<-(res$pvalue<0.05)*1
      net_1<-apply(net_1,c(1,2), function(x){ifelse(is.na(x),0,x)})
      lapl1<-as.dist(hamming.distance((net_1)))

      load(paste0("Results/",method,"_Day",df[j,1],".rdata"))
      net_2<-(res$pvalue<0.05)*1
      net_2<-apply(net_2,c(1,2), function(x){ifelse(is.na(x),0,x)})
      lapl2<-as.dist(hamming.distance((net_2)))

      mantel.pvalue_matrix[i,j]<- mantel.pvalue_matrix[j,i]<-mantel.rtest(lapl1,lapl2)$pvalue
    }
  }
  return(mantel.pvalue_matrix)
}

res<-mantel.pvalue(method="LUPINE", minday=2, maxday=10)

rownames(res)<-c(paste0("D",c(1,2,5, 6,7,9,12:14)))

combined_breaks <- c(seq(0, 0.05, 0.001), seq(0.0501,1,0.001))
combined_colors <- c(rev(colorRampPalette(brewer.pal(9,"RdPu")[1:5])(51)),
                     colorRampPalette(brewer.pal(9,"Blues"))(950))
combined_ramp <- colorRamp2(combined_breaks, combined_colors)
# Custom column annotation for label colors


col_ha <- HeatmapAnnotation(
  col_labels = anno_text(rownames(res),
                         rot = 360,
                         gp = gpar(col = c(rep('#009E73',3), rep('#F68B33',2), rep('#388ECC',4)),
                                   fontsize = 12)))


pdf("Figures/stattest_lupine.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(res, col=combined_ramp,
                        show_heatmap_legend=F,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F,
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep('#009E73',3), rep('#F68B33',2), rep('#388ECC',4))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", res[i, j]), x, y, gp = gpar(fontsize = 14))
                        })
par(op)
dev.off()
