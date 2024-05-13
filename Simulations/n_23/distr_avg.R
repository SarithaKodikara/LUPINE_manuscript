library(igraph)
library(mixOmics)
library(readr)
library(NetworkDistance)
library(rlist)
library(patchwork)
library(abind)

para<-readRDS("Sim_para/para.rds")
node_measures<-function(method="bPLS",iter,  day){
  if(!(day==1 &method=="bPLS")){
    load(paste0("Results/",method,"/iter",iter,"_d",day,".rdata"))
    
    if(method%in%c("SpiecEasi_mb","SpiecEasi_glasso")){
      net<-res$Graph
    }else{
      net<-(res$pvalue>0.95)*1
      net<-apply(net,c(1,2), function(x){ifelse(is.na(x),0,x)})
      
    }
    #net_ig<-graph_from_adjacency_matrix(net, mode = "undirected")
    return(net)
    
  }
  
}
True1<-para$binary.graph1_n
True2<-para$binary.graph1_h

dist_gdd<-function(method,iter){
  g_ls<-lapply(2:10,function(k){node_measures(method,iter,k)})
  
  g_ls_full<-list.append(g_ls,True1,True2)
  dst<-nd.gdd(g_ls_full, out.dist = FALSE)$D
  return(dst)
}

fit_df<-function(method, names, phase){
  dst_ls<-lapply(1:50, function(k){dist_gdd(method,k)})
  #dst<-Reduce("+",dst_ls)/50
  abind_ls <- do.call(abind, c(dst_ls, list(along=3)))
  dst<-apply(abind_ls, 1:2, mean)
  fit<-data.frame(cmdscale(dst, k=3))
  
  x<-fit[,1]
  y<-fit[,2]
  
  fit$name<-factor(names, level=names)
  fit$phase <-phase
  
  return(fit)
  
}

names1<-c(paste0("D", 2:10), "T[1]", "T[2]")
phase1<-factor(c(rep(c('Inferred Network 1','Inferred Network 2'), 
                     each=4),"Inferred Network 2",'True Network 1', 
                 'True Network 2'),
               levels = c('True Network 1', 'True Network 2',
                          'Inferred Network 1','Inferred Network 2'))


fit_bPLS<-fit_df("bPLS",names1, phase1)
fit_PCA<-fit_df("PCA",names1, phase1)
fit_sparcc<-fit_df("sparcc",names1, phase1)
fit_mb<-fit_df("SpiecEasi_mb",names1, phase1)
fit_glasso<-fit_df("SpiecEasi_glasso",names1, phase1)

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

# Shape
fit_bPLS$title <- "LUPINE"
fit_bPLS$name<-c(paste0("D[",2:10,"]"), paste0("T[",1:2,"]"))
(p1<-ggplot(fit_bPLS, aes(x=X1, y=X2,  color= phase)) + 
  geom_point(size=3, show.legend = FALSE)+ 
  labs(y = "MDS2", x = "MDS1")+
    theme_khakiwhite()+
    geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
    ylim(-1.5,1.5)+ xlim(-6,6)+facet_wrap(~title)+ 
  scale_colour_manual(values = c('#C2C2C2', '#009E73','#388ECC','#F68B33')))

fit_PCA$title <- "LUPINE_single"
fit_PCA$name<-c(paste0("D[",2:10,"]"), paste0("T[",1:2,"]"))
(p2<-ggplot(fit_PCA, aes(x=X1, y=X2, color= phase)) + 
  geom_point(size=3, show.legend = FALSE)+ 
  labs(y = "MDS2", x = "MDS1")+
    theme_khakiwhite()+
    geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
    ylim(-1.5,1.5)+ xlim(-6,6)+facet_wrap(~title)+ 
  scale_colour_manual(values = c('#C2C2C2', '#009E73','#388ECC','#F68B33')))

fit_sparcc$title = "sparcc"
fit_sparcc$name<-c(paste0("D[",2:10,"]"), paste0("T[",1:2,"]"))
(p3<-ggplot(fit_sparcc, aes(x=X1, y=X2, color= phase)) + 
  geom_point(size=3, show.legend = FALSE)+  
  labs(y = "MDS2", x = "MDS1")+
    theme_khakiwhite()+
    geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
    ylim(-1.5,1.5)+ xlim(-6,6)+facet_wrap(~title)+ 
  scale_colour_manual(values = c('#C2C2C2', '#009E73','#388ECC','#F68B33')))

fit_mb$title = "SpiecEasi_mb"
fit_mb$name<-c(paste0("D[",2:10,"]"), paste0("T[",1:2,"]"))
(p4<-ggplot(fit_mb, aes(x=X1, y=X2,  color= phase)) + 
  geom_point(size=3, show.legend = FALSE)+ 
  labs(y = "MDS2", x = "MDS1")+
    theme_khakiwhite()+
    geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
    ylim(-1.5,1.5)+ xlim(-6,6)+facet_wrap(~title)+ 
  scale_colour_manual(values = c('#C2C2C2', '#009E73','#388ECC','#F68B33')))

fit_glasso$title = "SpiecEasi_glasso"
fit_glasso$name<-c(paste0("D[",2:10,"]"), paste0("T[",1:2,"]"))
(p5<-ggplot(fit_glasso, aes(x=X1, y=X2, color= phase)) + 
  geom_point(size=3, show.legend = FALSE)+ 
  labs(y = "MDS2", x = "MDS1")+
    theme_khakiwhite()+
    geom_text(aes(label=name), parse=TRUE,hjust=0.5, vjust=1.5, show.legend = FALSE)+
  ylim(-1.5,1.5)+ xlim(-6,6)+facet_wrap(~title)+ 
  scale_colour_manual(values = c('#C2C2C2', '#009E73','#388ECC','#F68B33')))



pdf("Figures/n23_pcoa1_n.pdf",
    width=3, height=3)

p1

dev.off()


#legend

p1<-ggplot(fit_bPLS, aes(x=X1, y=X2, label = name, color= phase)) + 
  geom_point(size=3)+ 
  labs(y = "MDS2", x = "MDS1")+
  theme_khakiwhite()+ geom_text(hjust=0.5, vjust=1.5, show.legend = FALSE)+
  ylim(-1.5,1.5)+ xlim(-6,6)+facet_wrap(~title)+ 
  scale_colour_manual(values = c('#C2C2C2', '#009E73','#388ECC','#F68B33'))+
  theme(legend.direction="horizontal")
# Draw Only Legend without plot 
# Grab legend from gplot 
legend <- get_legend(p1)                     

# Create new plot window 
grid.newpage()                               

pdf("Figures/mds_legend.pdf",
    width=10, height=1)

grid.draw(legend)  

dev.off()

