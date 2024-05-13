library(ggplot2)
library(magrittr)
library(dplyr)
library(stringr)
library(tibble)
library(tidyr)

metafull_df <- read.csv("Data/Mu_etal_mouseVRE_MAPPING_FILE.txt",sep = "\t")

#filtering day 8 due to kissing OTU values for B2/B3
meta_DF<-data.frame(Sample.ID=paste0("X",metafull_df$X.SampleID),
                    Time=metafull_df$day_of_experiment,
                    Subject.ID=metafull_df$host_subject_id)%>%
  filter(., Subject.ID!="not applicable")%>%
  filter(.,Time!=8)

OTU_tsv <- read.csv("Data/Mu_etal_mouseVRE_FEATURE_TABLE_case_corrected.tsv",sep = "\t",skip = 1)

taxa_clasification<-read.csv( "Data/taxonomy.tsv",sep = "\t")
taxa_ordered<-taxa_clasification[order(taxa_clasification$Taxon),]

df_taxonomy<-data.frame(OTU=paste0("OTU_",1:3574),
                        cbind(Feature.ID=taxa_ordered$Feature.ID,
                              Taxon=taxa_ordered$Taxon,
                              str_split_fixed(taxa_ordered$Taxon, ";", 7)))%>%
              filter(V3=="k__Bacteria")
rownames(df_taxonomy)<-df_taxonomy$OTU

df_taxonomy$V4<-ifelse(df_taxonomy$V4=="", "p__",df_taxonomy$V4)
df_taxonomy$V5<-ifelse(df_taxonomy$V5=="", "c__",df_taxonomy$V5)
df_taxonomy$V6<-ifelse(df_taxonomy$V6=="", " o__",df_taxonomy$V6)
df_taxonomy$V7<-ifelse(df_taxonomy$V7=="", "f__",df_taxonomy$V7)
df_taxonomy$V8<-ifelse(df_taxonomy$V8=="", "g__",df_taxonomy$V8)
df_taxonomy$V9<-ifelse(df_taxonomy$V9=="", "s__",df_taxonomy$V9)



OTU <-OTU_tsv%>%slice(match(df_taxonomy$Feature.ID,X.OTU.ID))%>%
      left_join(., df_taxonomy[,1:2], by=join_by(X.OTU.ID==Feature.ID))
OTU_matrix<-t(OTU[,meta_DF$Sample.ID])
colnames(OTU_matrix)<-OTU[,110]

#checking if the two files are in same order
sum(rownames(OTU_matrix)!=meta_DF$Sample.ID)

#### filtering
all.fil.otus<-t(apply(OTU_matrix,1, function(x) x/sum(x))) %>%
  data.frame() %>%
  rownames_to_column('Sample.ID')%>%
  pivot_longer(cols = -Sample.ID,
               names_to = 'sOTUs',
               values_to = 'RelativeAbundance') %>% 
  mutate(RelativeAbundance = 100*RelativeAbundance) %>%
  left_join(meta_DF, by = 'Sample.ID')  %>% 
  group_by(sOTUs, Time) %>%
  summarise(MeanRA = mean(RelativeAbundance)) %>%
  filter(MeanRA >= 0.1) 


OTUtoKeep<- all.fil.otus%>%
  pull(sOTUs) %>%
  unique()  
#126 OTUs

taxanomy_filter<-df_taxonomy[OTUtoKeep,]
taxanomy_filter_ordered<-taxanomy_filter[order(taxanomy_filter$V6),]

OTUdata_array<-array(NA, c(9,126,10))
Lib_size<-array(NA, c(9,10))
colnames(OTUdata_array)<-taxanomy_filter_ordered$OTU
rownames(OTUdata_array)<-unique(meta_DF$Subject.ID)
index=0
for(i in unique(meta_DF$Time)){
  index=index+1
  tmp<-meta_DF[meta_DF$Time==i,]
  #Filtering out the OTU counts
  OTUdata_array[,,index]<-as.matrix(OTU_matrix[tmp$Sample.ID,taxanomy_filter_ordered$OTU])
  Lib_size[,index]<-rowSums(OTUdata_array[,,index])
  
}

##### Extract low abundance OTUs  across group and day

highabundance<-pivot_wider(all.fil.otus, names_from = sOTUs, values_from = MeanRA)%>%
  select(c("Time", taxanomy_filter_ordered$OTU))%>%
  arrange(Time)
OTU_l.abundance<-apply(X = is.na(highabundance),MARGIN = 1,
                       FUN = function(x) colnames(highabundance)[x])

save(OTUdata_array, file="data_filtered/OTUdata.rds")
save(Lib_size, file="data_filtered/Lib_size.rds")
save(taxanomy_filter_ordered, file="data_filtered/taxonomy.rds")
save(OTU_l.abundance, file="data_filtered/OTU_l.abundance.rds")
