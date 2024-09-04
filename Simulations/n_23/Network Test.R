library(graphsim)
library(ade4)
library(e1071)
library(abind)
library(RColorBrewer)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
set.seed(987)
mantel.pvalue<-function(method="bPLS",iter, minday=2, maxday=10){
  
  days= maxday-(minday-1)
  mantel.pvalue_matrix <- matrix(0, days, days)

  # Calculate pairwise mantel.pvalue
  for (i in 1:(days-1)) {
    for (j in (i+1):days) {
      
      day_vec<-minday:maxday
      
      load(paste0("Results/",method,"/iter",iter,"_d",day_vec[i],".rdata"))
      if(method%in%c("SpiecEasi_mb","SpiecEasi_glasso")){
        net_1<-res$Graph
      }else{
        net_1<-(res$pvalue>0.95)*1
      }
      net_1<-apply(net_1,c(1,2), function(x){ifelse(is.na(x),0,x)})
      lapl1<-as.dist(hamming.distance((net_1)))
      
      load(paste0("Results/",method,"/iter",iter,"_d",day_vec[j],".rdata"))
      if(method%in%c("SpiecEasi_mb","SpiecEasi_glasso")){
        net_2<-res$Graph
      }else{
        net_2<-(res$pvalue>0.95)*1
      }
      net_2<-apply(net_2,c(1,2), function(x){ifelse(is.na(x),0,x)})
      lapl2<-as.dist(hamming.distance((net_2)))
      if(sum(lapl1)!=0 & sum(lapl2)!=0 ){
        mantel.pvalue_matrix[i,j]<- mantel.pvalue_matrix[j,i]<-mantel.rtest(lapl1,lapl2)$pvalue
      }
      
    }
  }
    return(mantel.pvalue_matrix)
}


test_networks<-function(method){
  dst_ls<-lapply(1:50, function(k){mantel.pvalue(method,k)})
  abind_ls <- do.call(abind, c(dst_ls, list(along=3)))
  mean_pvalue <- apply(abind_ls, 1:2, mean)
  return(mean_pvalue)
}

pval_bPLS<-test_networks("bPLS")
pval_PCA<-test_networks("PCA")
pval_sparcc<-test_networks("SparCC")
pval_mb<-test_networks("SpiecEasi_mb")
pval_gl<-test_networks("SpiecEasi_glasso")
rownames(pval_PCA)<-paste0("D",2:10)
rownames(pval_bPLS)<-paste0("D",2:10)
rownames(pval_sparcc)<-paste0("D",2:10)
rownames(pval_mb)<-paste0("D",2:10)
rownames(pval_gl)<-paste0("D",2:10)


combined_breaks <- c(seq(0, 0.05, 0.001), seq(0.0501,1,0.001))
combined_colors <- c(rev(colorRampPalette(brewer.pal(9,"RdPu")[1:5])(51)),
                     colorRampPalette(brewer.pal(9,"Blues"))(950))
combined_ramp <- colorRamp2(combined_breaks, combined_colors)
# Custom column annotation for label colors

col_ha <- HeatmapAnnotation(
  col_labels = anno_text(paste0("D",2:10), 
                         rot = 360, 
                         gp = gpar(col = c(rep("yellow4",4), rep("lightsalmon4", 5)), 
                                   fontsize = 14)))
pdf("Figures/sim_n23_stattest_pca.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(pval_PCA, col=combined_ramp,
                        show_heatmap_legend=F,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F, 
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep("yellow4",4), rep("lightsalmon4", 5))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", pval_PCA[i, j]), x, y, gp = gpar(fontsize = 14))
                        })
par(op)
dev.off()

pdf("Figures/sim_n23_stattest_bPLS.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(pval_bPLS, col=combined_ramp,
                        show_heatmap_legend=T,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F, 
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep("yellow4",4), rep("lightsalmon4", 5))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", pval_bPLS[i, j]), x, y, gp = gpar(fontsize = 14))
                        })
par(op)
dev.off()


pdf("Figures/sim_n23_stattest_sparcc.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(pval_sparcc, col=combined_ramp,
                        show_heatmap_legend=F,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F, 
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep("yellow4",4), rep("lightsalmon4", 5))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", pval_sparcc[i, j]), x, y, gp = gpar(fontsize = 14))
                        })
par(op)
dev.off()


pdf("Figures/sim_n23_stattest_mb.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(pval_mb, col=combined_ramp,
                        show_heatmap_legend=F,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F, 
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep("yellow4",4), rep("lightsalmon4", 5))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", pval_mb[i, j]), x, y, gp = gpar(fontsize = 14))
                        })
par(op)
dev.off()


pdf("Figures/sim_n23_stattest_gl.pdf",
    width=6, height=6)
op <- par(mar = rep(0, 4))
ComplexHeatmap::Heatmap(pval_gl, col=combined_ramp,
                        show_heatmap_legend=F,
                        border=0,  name = "p-value",
                        cluster_rows = F, cluster_columns = F, 
                        top_annotation = col_ha,
                        row_names_gp = gpar(col = c(rep("yellow4",4), rep("lightsalmon4", 5))),
                        cell_fun = function(j, i, x, y, width, height, fill) {
                          grid.text(sprintf("%.2f", pval_gl[i, j]), x, y, gp = gpar(fontsize = 14))
                        })
par(op)
dev.off()
