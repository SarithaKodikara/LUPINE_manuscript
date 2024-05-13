
library(dplyr) #For '%>%'
library(tibble) # For 'rownames_to_column'
library(tidyr) # For 'pivot_longer'
library(stringr) # For 'str_split_fixed'

#### Reading Data ####
load('Data_org/BEME.RData') 
#sum(data.raw==0)

####Removing the offset added to the data
data.raw<-data.raw-1
sum(data.raw==0)


###filtering Biological samples, very high library size [sort(rowSums(data.raw))] and the initial timepoint
#filtered_bio_OTU<-data.raw%>%  data.frame()%>%filter(sample.info$Rep=="Biological" )
#df_libsize<-data.frame(rowNames = rownames(filtered_bio_OTU), lib.size=rowSums(filtered_bio_OTU))
#hist(df_libsize$lib.size, xlab = "Library Size", main = "")

filtered_bio_OTU<-data.raw%>%  data.frame()%>%
  filter(sample.info$Rep=="Biological" & rownames(.)!="83M" &sample.info$Day!="A")

df_libsize<-data.frame(sample_name = rownames(filtered_bio_OTU), lib.size=rowSums(filtered_bio_OTU))

#### filtering data by minimum reads ans minimum occurance(i.e., % of zeros) ####
sample.info<- sample.info%>%rownames_to_column('sample_name')
OTUtoKeep<-t(apply(filtered_bio_OTU,1, function(x) x/sum(x))) %>%
            data.frame() %>%
                      rownames_to_column('sample_name')%>%
                      pivot_longer(cols = -sample_name,
                                   names_to = 'sOTUs',
                                   values_to = 'RelativeAbundance') %>% 
            mutate(RelativeAbundance = 100*RelativeAbundance) %>%
            left_join(sample.info, by = 'sample_name')  %>% 
            group_by(sOTUs, Diet, Day) %>%
            summarise(MeanRA = mean(RelativeAbundance)) %>%
            filter(MeanRA >= 0.1) %>%
            pull(sOTUs) %>%
            unique()  
#212 OTUs

filtered_OTU<-filtered_bio_OTU[,OTUtoKeep]

### filtering meta information ####
filtered_sample.info<-sample.info[sample.info$sample_name %in%
                                    rownames(filtered_OTU),] %>%left_join(.,df_libsize)


### filter taxanomy information ####

df_taxonomy<-data.frame(cbind(colnames(data.raw),str_split_fixed(taxonomy, ";", 7)))
df_taxonomy$X5<-ifelse(df_taxonomy$X5=="", " o__",df_taxonomy$X5)
df_taxonomy$X6<-ifelse(df_taxonomy$X6=="", " f__",df_taxonomy$X6)
filtered_taxonomy<-df_taxonomy[df_taxonomy$X1%in%colnames(filtered_OTU),]
rownames(filtered_taxonomy)<-filtered_taxonomy$X1
filtered_taxonomy$X5<-factor(filtered_taxonomy$X5)


#Finding mouse that have a complete data across time

mouseID_norm<-sort(Reduce(intersect, 
       list(filtered_sample.info$MouseID[filtered_sample.info$Diet=="Normal" & filtered_sample.info$Day=="Day0"],
            filtered_sample.info$MouseID[filtered_sample.info$Diet=="Normal" & filtered_sample.info$Day=="Day1"],
            filtered_sample.info$MouseID[filtered_sample.info$Diet=="Normal" & filtered_sample.info$Day=="Day4"],
            filtered_sample.info$MouseID[filtered_sample.info$Diet=="Normal" & filtered_sample.info$Day=="Day7"])))  

mouseID_hfhs<-sort(Reduce(intersect, 
                          list(filtered_sample.info$MouseID[filtered_sample.info$Diet=="HFHS" & filtered_sample.info$Day=="Day0"],
                               filtered_sample.info$MouseID[filtered_sample.info$Diet=="HFHS" & filtered_sample.info$Day=="Day1"],
                               filtered_sample.info$MouseID[filtered_sample.info$Diet=="HFHS" & filtered_sample.info$Day=="Day4"],
                               filtered_sample.info$MouseID[filtered_sample.info$Diet=="HFHS" & filtered_sample.info$Day=="Day7"])))  

#Sorting OTUs based on total counts across order
OTU_healthy<-data.frame(OTU=colnames(filtered_OTU),
                        Sum=colSums(filtered_OTU[filtered_sample.info$sample_name[filtered_sample.info$Diet=="Normal"],]),
                        Order=filtered_taxonomy[colnames(filtered_OTU),5])
OTU_healthy_sort<-arrange(OTU_healthy,Order,desc(Sum))


OTUdata_Normal<-array(NA, c(23,212,4))
colnames(OTUdata_Normal)<-OTU_healthy_sort$OTU
rownames(OTUdata_Normal)<-paste0("M_",mouseID_norm)
  
OTUdata_HFHS<-array(NA, c(23,212,4))
colnames(OTUdata_HFHS)<-OTU_healthy_sort$OTU
rownames(OTUdata_HFHS)<-paste0("M_",mouseID_hfhs)

Lib_Normal<-matrix(NA, nrow = 23, ncol=4)
Lib_HFHS<-matrix(NA, nrow = 23, ncol=4)

for(i in 1:4){
  
  selectD<-ifelse(i==1,"Day0", ifelse(i==2, "Day1", ifelse(i==3, "Day4", "Day7")))
  
  #Ordering mouse across all timepoints
  healthy_tmp<-filtered_sample.info[filtered_sample.info$MouseID%in%mouseID_norm & 
                         filtered_sample.info$Day==selectD,]
  healthy_tmp.ord<-healthy_tmp[order(healthy_tmp$MouseID),]
  #Filtering out the OTU counts
  OTUdata_Normal[,,i]<-as.matrix(filtered_OTU[healthy_tmp.ord$sample_name,OTU_healthy_sort$OTU])
  Lib_Normal[,i]<-healthy_tmp.ord$lib.size
  
  hfhs_tmp<-filtered_sample.info[filtered_sample.info$MouseID%in%mouseID_hfhs & 
                                      filtered_sample.info$Day==selectD,]
  hfhs_tmp.ord<-hfhs_tmp[order(hfhs_tmp$MouseID),]
  OTUdata_HFHS[,,i]<-as.matrix(filtered_OTU[hfhs_tmp.ord$sample_name,OTU_healthy_sort$OTU])
  Lib_HFHS[,i]<-hfhs_tmp.ord$lib.size
}
  
save(OTUdata_HFHS, file="Data_filtered/OTUdata_HFHS.rds")
save(OTUdata_Normal, file="Data_filtered/OTUdata_Normal.rds")
save(Lib_HFHS, file="Data_filtered/Lib_HFHS.rds")
save(Lib_Normal, file="Data_filtered/Lib_Normal.rds")
save(filtered_sample.info, file="Data_filtered/filtered_sample.info.rds")
save(filtered_taxonomy, file="Data_filtered/filtered_taxonomy.rds")

