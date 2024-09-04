
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
library(abind)
library(tidyr)
library(dplyr)

Plot_pvalues<- function(median_mt, lower_mt,upper_mt, taxa_colors, topNumber= 10, taxanames){
  #Extract the upper triangle values (excluding the diagonal)
  upper_triangle_indices <- which(upper.tri(median_mt, diag = FALSE), arr.ind = TRUE)
  estimate_upper_triangle <- median_mt[upper_triangle_indices]


  # Find the indices of the 10 smallest upper triangle values
  smallest_indices <- order(estimate_upper_triangle, decreasing = FALSE)[1:topNumber]

  # Extract the corresponding values from lower and upper matrices
  selected_indices <- upper_triangle_indices[smallest_indices, ]
  lower_vals <- lower_mt[selected_indices]
  upper_vals <- upper_mt[selected_indices]
  estimate_vals <- median_mt[selected_indices]



  #  Create a data frame for easier handling
  data <- data.frame(
    index =1:topNumber,
    estimate = estimate_vals,
    lower = lower_vals,
    upper = upper_vals,
    row = selected_indices[, 1],      # Row indices
    col = selected_indices[, 2],      # Column indices
    row_col =paste0(taxanames[selected_indices[, 1]],"--",taxanames[selected_indices[, 2]])
    #row_col = paste0( "Taxa",selected_indices[, 1], " -- Taxa", selected_indices[, 2])
  )


  # Convert row indices to factor for plotting
  data$index <- factor(data$index)


  data$row_color <-taxa_colors[data$col]
  data$col_color <-taxa_colors[data$row]

  # Step 6: Plot using ggplot2 with flipped axes and remove boxes
  p<-ggplot(data, aes(x = index, y = log(estimate))) +
    geom_errorbar(aes(ymin = log(lower), ymax = log(estimate)), width = 0.2, size=1.5, col= data$col_color) +
    geom_errorbar(aes(ymin = log(estimate), ymax = log(upper)), width = 0.2, size=1.5, col= data$row_color) +
    geom_point(size = 3) +
    coord_flip() + # Flip the axes
    xlab("") +
    ylab("log(pvalue)") +
    theme_minimal() +
    theme(
      panel.border = element_blank(),     # Remove panel border
      axis.line = element_blank(),        # Remove axis lines
      axis.ticks = element_blank(),       # Remove axis ticks
      plot.title = element_text(hjust = 0.5), # Center the title
      panel.grid.major.x = element_blank() ,
      panel.grid.minor.x = element_blank(),
      legend.position = 'none',
      axis.text.y = element_text(size = 12)
    )+ylim(-80,0)+
    scale_x_discrete(labels=as.character(data$row_col))+
    geom_hline(yintercept = log(0.05), color = "grey66", linetype = "dashed", size = 1)

  return(p)
}



load("Results/LUPINE_boot.rdata")
load("data_modified/taxonomy.rds")
load("data_modified/data_PPT.rds")
Day0<-data_PPT[,,1]
taxa_info<-taxonomy_order[match(colnames(Day0), taxonomy_order$Taxa),]

taxa_genus<-gsub("_unclassified","*",gsub("g__","",taxa_info$X6))

taxa_colors <- rep(c( "yellow", "green",
     "darkgreen",  "pink", "hotpink",  "royalblue",
     "darkred", "firebrick2",  "bisque4", "orange" ,
     "cyan3", "steelblue4",  "darkolivegreen3", "plum1",
     "coral2", "purple",  "aquamarine3", "khaki4",
     "indianred", "mediumvioletred",  "tan4",
     "turquoise4",  "wheat", "springgreen", "black"),
  c(summary(factor(taxa_info$X6))))

##PPT
pvalues<-inf_upper[[1]] %>%
  abind(., along=3)

median_mt<-apply(pvalues, c(1,2), median)
lower_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.025))
upper_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.975))


pdf("Figures/Prediabetic_CI_PPT.pdf",
    width=6, height=3)
op <- par(mar = rep(0, 4))
Plot_pvalues(median_mt,lower_mt, upper_mt,taxa_colors,15, taxa_genus)
par(op)
dev.off()


##MDE
pvalues<-inf_upper[[2]] %>%
  abind(., along=3)

median_mt<-apply(pvalues, c(1,2), median)
lower_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.025))
upper_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.975))


pdf("Figures/Prediabetic_CI_MDE.pdf",
    width=6, height=3)
op <- par(mar = rep(0, 4))
Plot_pvalues(median_mt,lower_mt, upper_mt,taxa_colors,15, taxa_genus)
par(op)
dev.off()
