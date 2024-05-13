#https://github.com/gerberlab/mitre/tree/master/mitre/example_data/david

library(tidyverse)
library(compositions)
library(splines)
library(ggplot2)
library(dplyr)
library(magrittr)
library(phyloseq)
library(fs)
library(glue)
library(seqinr)
library(dada2)


reads <- read_csv( "data/abundance.csv") |>
  dplyr::rename(sample = `...1`) |>
  column_to_rownames("sample") |>
  as.matrix()

samples <- read_csv("data/sample_metadata.csv", col_names = c("sample", "subject", "time")) |>
  mutate(group = str_extract(subject, "[A-z]+"))

## Taxonomy
#taxonomy_full <- assignTaxonomy("data/sequence_key.fa", "data/silva_nr99_v138.1_wSpecies_train_set.fa.gz")
#taxonomy_full<-cbind(OTU=colnames(reads),taxonomy_full)
#save(taxonomy_full, file="data_modified/taxonomy_full.rds")
load("data_modified/taxonomy_full.rds")

all.fil.otus<-t(apply(reads,1, function(x) x/sum(x))) %>%
  data.frame() %>%
  rownames_to_column('sample')%>%
  pivot_longer(cols = -sample,
               names_to = 'sOTUs',
               values_to = 'RelativeAbundance') %>% 
  mutate(RelativeAbundance = 100*RelativeAbundance) %>%
  left_join(samples, by = 'sample')  %>% 
  group_by(sOTUs, group, time) %>%
  summarise(MeanRA = mean(RelativeAbundance)) %>%
  filter(MeanRA >= 0.1) 


OTUtoKeep<- all.fil.otus %>%
  pull(sOTUs) %>%
  unique()  
#311 OTUs

### filter taxanomy information ####
taxonomy_fil<-taxonomy_full%>% data.frame%>% filter(OTU%in% OTUtoKeep)
taxonomy_fil$Order<-ifelse(is.na(taxonomy_fil$Order), "o__", paste0("o__", taxonomy_fil$Order))
taxanomy_filter_ordered<-taxonomy_fil[order(taxonomy_fil$Order),]


reads_fil<-reads[,taxanomy_filter_ordered$OTU]
reads_norm<-apply(reads_fil,1,compositions::clr)%>%t() |>
  data.frame() |>
  rownames_to_column("sample")

combined_all<-samples%>%
  left_join(reads_norm , by = 'sample') %>%
  mutate_at(vars(subject, time), factor)

#NA values for missing days
#combined_full<-merge(combined_all, expand.grid(lapply(combined_all[c("subject","time")], levels)), all.y=TRUE)


ggplot(combined_all, aes(x=time, y=subject, group=group, color=group)) +
  geom_point() 

#### Interpolate missing timepoints 
interpolate_fun<-function(data, subject){
  
  data_fil<-data[data$subject==subject,]%>%
    arrange(time)%>%data.frame()
  #Converting factor into integer
  time_i<-as.integer(data_fil$time)
  data_i<-data_fil[, -c(1:4)]
  
  new_data<-sapply(1:311, function(i){
    full_y<-rep(NA, 16)
    all_time  <- 1:16
    y_interpolate <- spline(x=time_i,  y=data_i[,i], xout = all_time)$y      
    full_y[all_time]<-y_interpolate
    full_y
    #plot(time_i,  data_i[,i])
    #lines(predict(sspline, all_times), col = 3)
  })%>%
    set_colnames(taxanomy_filter_ordered$OTU)
  new_df<-data.frame(subject=rep(data_fil$subject[1],16),
                     group=rep(data_fil$group[1], 16),
                     time=-5:10)%>%
    cbind(new_data)
  
  return(new_df)
  
}

#interpolate all subjects |> convert to data frame |> extract days from -4 to 10
reads_inter<-sapply(levels(combined_all$subject),
                    function(j){interpolate_fun(combined_all,j)}, simplify = FALSE)%>%
  plyr::ldply( data.frame)%>%
  filter(time %in% c(-4:10))


##### Extract low abundance OTUs  across group and day
highabundance<-pivot_wider(all.fil.otus, names_from = sOTUs, values_from = MeanRA)%>%
  dplyr::select(c("group", "time", all_of(taxanomy_filter_ordered$OTU)))%>%
  filter(time %in% c(-4:10))

OTU_l.abundance<-apply(X = is.na(highabundance),MARGIN = 1,
                       FUN = function(x) colnames(highabundance)[x])


low_animal<-OTU_l.abundance[1:15]
low_plant<- OTU_l.abundance[16:30] 

# create matrix for each group

OTUdata_plant<-array(NA, c(10,311,15))
colnames(OTUdata_plant)<-taxanomy_filter_ordered$OTU
rownames(OTUdata_plant)<-c("Plant1","Plant2","Plant3","Plant4","Plant5",
                           "Plant6","Plant7","Plant8","Plant9","Plant10")

OTUdata_animal<-array(NA, c(10,311,15))
colnames(OTUdata_animal)<-taxanomy_filter_ordered$OTU
rownames(OTUdata_animal)<-c("Animal1","Animal2","Animal3","Animal4","Animal5",
                            "Animal6","Animal8","Animal9","Animal10","Animal11")



for(i in -4:10){
  index=i+5
  plat_tmp<-reads_inter[reads_inter$time==i &
                          reads_inter$group=="Plant",]%>%
    arrange(subject)
  OTUdata_plant[,,index]<-as.matrix(plat_tmp[,taxanomy_filter_ordered$OTU])
  
  animal_tmp<-reads_inter[reads_inter$time==i &
                            reads_inter$group=="Animal",]%>%
    arrange(subject)
  OTUdata_animal[,,index]<-as.matrix(animal_tmp[,taxanomy_filter_ordered$OTU])
}

filtered_sample.info<-reads_inter[,1:4]



save(taxanomy_filter_ordered, file="data_modified/taxonomy_fil.rds")
save(OTUdata_plant, file="data_modified/OTUdata_plant.rds")
save(OTUdata_animal, file="data_modified/OTUdata_animal.rds")
save(filtered_sample.info, file="data_modified/filtered_sample.info.rds")
save(low_animal, file="data_modified/low_animal.rds")
save(low_plant, file="data_modified/low_plant.rds")

