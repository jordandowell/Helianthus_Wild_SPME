# data processing
#libraries used
library("readxl")

#i need to incorporate batch folder reading 

#create empty  lists
#create empty dataframe to store batch and sample information 
SampleIDtable<-data.frame()


#start for loop K for batchfiles 

#batch files
BATCHFILES<-list.files("Data/Jessica_SPME/")

for (k in 1:length(BATCHFILES)) {
  


#import sample mass and metadata

#this is for just batch one
#get file name list
Shimadzufiles <-
  list.files(path =paste0("Data/Jessica_SPME/",BATCHFILES[k],"/"), pattern = "ASCII")


#for loop paste ID and batch identifier 
i<-0
for (i in 1:length(Shimadzufiles)) {

primarytable<-read.table(
  file = paste0("Data/Jessica_SPME/",BATCHFILES[k],"/", Shimadzufiles[i]),
  skip = 17,
  nrows = 2
)

primaryinfo<-c(paste0("Data/Jessica_SPME/",BATCHFILES[k],"/", Shimadzufiles[i]),BATCHFILES[k],Shimadzufiles[i],primarytable[1,3],primarytable[2,3])
SampleIDtable<-rbind(SampleIDtable,primaryinfo)
}

}




#add proper column names to make things easier
colnames(SampleIDtable)<-c("filepath","batchname","SHIMADZUFILENAME","BATCHID","SAMPLEID")

k<-0
BATCHMETA<-data.frame()
for (k in 1:length(BATCHFILES)) {
#Sheet needs to be batch file number 
Batch1_Meta <-
  read_excel("Data/Jessica SPME 1_24_2021.xlsx", sheet = k)
BATCHMETA<-rbind(BATCHMETA,Batch1_Meta)
}

#bind file name to sample mass and meta data

Batch_Meta_SampleID<-(cbind(SampleIDtable,BATCHMETA))


#remove blanks from file list
Batch_Meta_SampleID_cleaned<-Batch_Meta_SampleID[!Batch_Meta_SampleID$Species=='blank',]


View(Batch_Meta_SampleID_cleaned$filepath[2])

length(unique(Batch_Meta_SampleID_cleaned$filepath))

#we want to read in files with the SampleID matching matching the files we want
#we are analysing compounds with at least a 70% similiarity score to NIST.
#this ensure we are assessing common compounds and how they are shifting.
# we are also pulling out unidentified peaks so there will be multiple data sets
#for loop l over shimadzu files
#start here

#get the list of the number of samples 
NumberOfSamples<-length(unique(Batch_Meta_SampleID_cleaned$filepath))


l<-0
k<-0
#create emptylist for the for loops
RT.Peakarea_list<-list()
RT.Peakarea_list_known<-list()
Peak_labels<-list()

for (k in 1:NumberOfSamples) {
  print(k)
 # Shimadzufiles <-
#    list.files(path =paste0("Data/Jessica_SPME/",BATCHFILES[k],"/"), pattern = "ASCII")
  
#for(l in 1:length(Shimadzufiles)){
#print(Shimadzufiles[l])    


#read in compound search results

#read in full data file
Trial<-readLines(Batch_Meta_SampleID_cleaned$filepath[k])
#find line where similarity starts & ends
table.Start<-which(Trial=="[MS Similarity Search Results for Spectrum Process Table]",arr.ind = TRUE)
table.End<-which(Trial=="[MS Chromatogram]",arr.ind = TRUE)


#save as data frame

datastream <-Trial[(table.Start+3):(table.End-1)]
SimilaritySearch_1<-(read.table(text=datastream,fill=T,sep = "\t",quote = ""))

(datastream[6])
#find the start of the MC peak table and save at dataframe
MC.table.Start<-which(Trial=="[MC Peak Table]",arr.ind = TRUE)
MC.table.End<-which(Trial=="[Spectrum Process Table]",arr.ind = TRUE)
datastream2<-Trial[(MC.table.Start+4):(MC.table.End-1)]
MC.Peak.Table<-(read.table(text=datastream2,fill=T,sep = "\t",quote = ""))
#create a table with just retention time and peak area
MC.RT.Peakarea<-MC.Peak.Table[,c(2,6)]


#pull out only compounds with known identity

MC.Peak.Table_Known<-MC.Peak.Table[MC.Peak.Table$V1%in%unique(SimilaritySearch_1$V1),]
#create a table with just retention time and peak area
MC.RT.Peakarea_Known<-MC.Peak.Table_Known[,c(2,6)]


#table of labeled peak number
Peak_Label_number<- as.data.frame(unique(SimilaritySearch_1$V1))
colnames(Peak_Label_number)<- paste0(Batch_Meta_SampleID_cleaned$filepath[k])



#add to lists for labels, rt for known and unknown peaks
RT.Peakarea_list[[length(RT.Peakarea_list)+1]]<-MC.RT.Peakarea
RT.Peakarea_list_known[[length(RT.Peakarea_list_known)+1]]<-MC.RT.Peakarea_Known
Peak_labels[length(Peak_labels)+1]<-Peak_Label_number


#for loop should end here 

}




#create dataframes that allow for unequal sizes for each data frame padd the end with 0's 


#cbind fill function from (https://gist.github.com/abelsonlive/4112423)
cbind.fill<-function(...){
  nm <- list(...) 
  nm<-lapply(nm, as.matrix)
  n <- max(sapply(nm, nrow)) 
  do.call(cbind, lapply(nm, function (x) 
    rbind(x, matrix(, n-nrow(x), ncol(x))))) 
}

#peaks with potential IDs
RT.Peakarea_list_known_Fulldata<-data.frame()
for (i in 1:length(RT.Peakarea_list_known)) {
  RT.Peakarea_list_known_Fulldata<-cbind.fill(RT.Peakarea_list_known_Fulldata,RT.Peakarea_list_known[[i]])
}

#all peaks

RT.Peakarea_list_Fulldata<-data.frame()
for (i in 1:length(RT.Peakarea_list)) {
  RT.Peakarea_list_Fulldata<-cbind.fill(RT.Peakarea_list_Fulldata,RT.Peakarea_list[[i]])
}


#peak label list

Peak_labels_list_Fulldata<-data.frame()
for (i in 1:length(Peak_labels)) {
  Peak_labels_list_Fulldata<-cbind.fill(Peak_labels_list_Fulldata,Peak_labels[[i]])
}


#save each dataframe to the output and make files necessary for GCalignR in Excel 

getwd()


write.csv(RT.Peakarea_list_known_Fulldata,"Output/OutputData/RT.Peakarea_KnownCompounds.csv")
write.csv(RT.Peakarea_list_Fulldata,"Output/OutputData/RT.Peakarea_AllCompounds.csv")
write.csv(Peak_labels_list_Fulldata,"Output/OutputData/PeakLabels_KnownCompounds.csv")
#write modified metadata
write.csv(Batch_Meta_SampleID_cleaned,"Output/OutputData/Sample_Metadata.csv")


