library(dplyr)
library(tidyr)




Peak_labels_annotation_list_Fulldata<-data.frame()
for (i in 1:length(PeaksWithLabes_list)) {
  
  PeaksWithLabes_list[[i]][,8]<-Batch_Meta_SampleID_cleaned[i,1]
  Peak_labels_annotation_list_Fulldata<-rbind(Peak_labels_annotation_list_Fulldata,PeaksWithLabes_list[[i]])
 # Peak_labels_annotation_list_Fulldata<-cbind.fill(Peak_labels_annotation_list_Fulldata,Batch_Meta_SampleID_cleaned[i,1])
}
dim(Peak_labels_annotation_list_Fulldata)

View(Peak_labels_annotation_list_Fulldata)

View(PeaksWithLabes_list[[1]])

#start of for loop 

#get list of unique compounds
listofallPotentialCMPDS<-rbind(Peak_labels_annotation_list_Fulldata[,seq(from=7,to=1728,8)])


Unique_listofallPotentialCMPDS<-as.data.frame(unique(as.vector(as.matrix(listofallPotentialCMPDS))))
colnames(Unique_listofallPotentialCMPDS)<-"COMPOUND"


#make a list of metadata names
Batch_Meta_SampleID_cleaned$SAMPLEMeta <- paste(Batch_Meta_SampleID_cleaned$Species,Batch_Meta_SampleID_cleaned$Treatment,Batch_Meta_SampleID_cleaned$Replicate,Batch_Meta_SampleID_cleaned$Mass,sep = "_")



for (i in 1:length(PeaksWithLabes_list)) {
 
  Peakinfo<- PeaksWithLabes_list[[i]][,c(2,3,7)]
  df<- Peakinfo%>%
    group_by(V5) %>%
    mutate(SAMPLERT = paste0(Batch_Meta_SampleID_cleaned$SAMPLEMeta[i],'_RT', row_number())) %>%
    mutate(SAMPLEArea = paste0(Batch_Meta_SampleID_cleaned$SAMPLEMeta[i],'_Area', row_number()))%>%
    spread(SAMPLERT, V2.x) %>%
    spread(SAMPLEArea,V6.x)
  
  df<-df %>% group_by(V5) %>% summarise_each(funs(sum(., na.rm = TRUE))) 
  Unique_listofallPotentialCMPDS<-merge(Unique_listofallPotentialCMPDS,df, by.x = "COMPOUND",by.y = "V5",all = T)
  
  
   
}


write.csv(Unique_listofallPotentialCMPDS,"Output/OutputData/AlignedbyIdentity.csv")


#i removed all silayted products and nonsense molecules 

FullCompounds<- read.csv("Output/OutputData/AlignedbyIdentity_NAs.csv")


View(FullCompounds)


#subset based on RT time
df1 <- FullCompounds[ , grepl( "RT" , names( FullCompounds ) ) ]
df2<-cbind(FullCompounds[,1],df1)
#change 0 to NA to get average RT
df2[df2 ==0]<-NA
mean_values <- cbind(df2[,1],rowMeans(df2[,2:ncol(df2)], na.rm = TRUE))

View(df2)

TDF<-t(df2)
colnames(TDF)<-TDF[1,]
TDF<-TDF[-1,]
TDF<-data.frame(apply(TDF, 2, as.numeric))

library(psych)
SummaryofCompounds<-describe(TDF,na.rm = TRUE)

write.csv(cbind(df2[,1],SummaryofCompounds),"Output/OutputData/RTsummaries.csv")
View(SummaryofCompounds)

View(TDF)

length(unique(FullCompounds[,1])) == nrow(FullCompounds)

which(duplicated(FullCompounds[,1]))


