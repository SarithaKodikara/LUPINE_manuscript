library(readxl)
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

raw_data<- read_excel("data/41467_2023_41042_MOESM9_ESM.xlsx", sheet = "gut species")

raw_data<- raw_data %>% 
              mutate(., SampleID=paste0("sample", 1:dim(raw_data)[1]), .before=1) %>% 
              mutate(., PatientID=paste0("patient_", `Participant ID`), .before=2) %>%
              mutate(., Time= if_else(`Time Point`=="Pre-intervention",0,1), .before=3)

sample.info<- raw_data[,1:4]

bac_names<- grep("k__Bacteria", colnames(raw_data), value = TRUE, ignore.case = FALSE)

bac_data<-raw_data[, bac_names]  

filtered_bac_data<- bac_data[, colMeans(!is.na(bac_data))>0.75]
filtered_bac_data[is.na(filtered_bac_data)]<-0

taxonomy<- data.frame(Feature.ID=colnames(filtered_bac_data),
                      str_split_fixed(colnames(filtered_bac_data), "\\|", 8))
taxonomy_order<-cbind(Taxa=paste0("Taxa", 1:dim(taxonomy)[1]),taxonomy[order(taxonomy$X6),])


filtered_bac_ordered<-filtered_bac_data[,taxonomy_order$Feature.ID]
colnames(filtered_bac_ordered)<-taxonomy_order$Taxa


data_MDE<-array(NA, c(96,91,2))
colnames(data_MDE)<-taxonomy_order$Taxa
rownames(data_MDE)<- sample.info%>% filter(.,Diet=="MED diet") %>% dplyr::select(PatientID)%>% unique()%>% pull()

data_PPT<-array(NA, c(93,91,2))
colnames(data_PPT)<-taxonomy_order$Taxa
rownames(data_PPT)<- sample.info%>% filter(.,Diet=="PPT diet") %>% dplyr::select(PatientID)%>% unique()%>% pull()

for(i in 0:1){
  index=i+1
  day=i
  data_MDE[,,index]<-as.matrix(filtered_bac_ordered[sample.info$Time==day &
                                         sample.info$Diet=="MED diet",])
  data_PPT[,,index]<-as.matrix(filtered_bac_ordered[sample.info$Time==day &
                                         sample.info$Diet=="PPT diet",])
}

save(taxonomy_order, file="data_modified/taxonomy.rds")
save( data_MDE, file="data_modified/data_MDE.rds")
save(data_PPT, file="data_modified/data_PPT.rds")
save(sample.info, file="data_modified/sample.info.rds")
