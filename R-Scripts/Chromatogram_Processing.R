# data processing
#libraries used
library("readxl")

#create dataframes to store full metadata and chromatograms

Full_Chromatogram_sets <- data.frame(row.names = 1:2840)
Full_Batch_Metadata <- data.frame()


#get number of batch folders
batchFolders <- list.files("Data/Jessica_SPME", full.names = T)

for (j in 1:length(batchFolders)) {
  #import sample mass and metadata
  Batch1_Meta <-
    read_excel("Data/Jessica SPME 1_24_2021.xlsx", sheet = j)
  #overwrite and keep only the first 5 columns of information
  Batch1_Meta <- Batch1_Meta[, 1:5]
  #this is for just batch one
  #get file name list
  Shimadzufiles <-
    list.files(path = batchFolders[j], pattern = "ASCII")
  #create empty dataframe to store batch and sample information
  SampleIDtable <- data.frame()
  
  
  #for loop paste ID and batch identifier
  i <- 0
  for (i in 1:length(Shimadzufiles)) {
    primarytable <- read.table(
      file = paste0(batchFolders[j], "/", Shimadzufiles[i]),
      skip = 17,
      nrows = 2
    )
    
    primaryinfo <-
      c(Shimadzufiles[i], primarytable[1, 3], primarytable[2, 3])
    SampleIDtable <- rbind(SampleIDtable, primaryinfo)
  }
  #order data frame so that it's in the same order as the batch meta data
  
  SampleIDtable <- SampleIDtable[order(SampleIDtable[, 3]), ]
  #add proper column names to make things easier
  colnames(SampleIDtable) <- c("SHIMADZUFILENAME", "BATCHID", "SAMPLEID")
  
  
  #bind file name to sample mass and meta data
  
  Batch1_Meta_SampleID <- (cbind(SampleIDtable, Batch1_Meta))
  
  
  #remove blanks from file list
  Batch1_Meta_SampleID_cleaned <-
    Batch1_Meta_SampleID[!Batch1_Meta_SampleID$Species == 'blank', ]
  
  #we want to read in files with the SampleID matching matching the files we want
  #we are analysing the chromatograms of samples
  
  
  
  
  #create an empty data frame
  Chromatogram_sets <- data.frame(row.names = 1:2840)
  #read in compound search results
  for (i in 1:dim(Batch1_Meta_SampleID_cleaned)[1]) {
    #read in full data file
    Trial <-
      readLines(paste0(
        batchFolders[j],
        "/",
        Batch1_Meta_SampleID_cleaned$SHIMADZUFILENAME[i]
      ))
    #find line where chromatogram starts & ends
    table.Start <- which(Trial == "[MS Chromatogram]", arr.ind = TRUE)
    #the table ends before the appearance of the first MS Spectrum
    table.End <- which(Trial == "[MS Spectrum]", arr.ind = TRUE)[1]
    #create a column name indicating batch# and file name
    Column_Name <-
      paste0("Batch_",
             j,
             "_",
             Batch1_Meta_SampleID_cleaned$SHIMADZUFILENAME[i])
    #save as data frame
    datastream <- Trial[(table.Start + 7):(table.End - 1)]
    Chromatogram_1 <- (read.table(
      text = datastream,
      fill = T,
      sep = "\t"
    ))
    #trim to remove scans before 5.8
    #* one sample started the MS @5 minute but there was nothing until about six minutes
    Chromatogram_2 <- Chromatogram_1[!(Chromatogram_1$V1 < 5.8), ]
    
    print(j)
    
    Chromatogram_sets[, Column_Name] <- Chromatogram_2[, 2]
    
    #add batch and file as columnname
    colnames(Chromatogram_sets)[i] <- Column_Name
    
    
  }
  #add time as rownames
  rownames(Chromatogram_sets) <- Chromatogram_2[, 1]
  
  
  
  
  
  Batch1_Meta_SampleID_cleaned$SampleBatchID <-
    colnames(Chromatogram_sets)
  
  
  #bind batch data to the full dataset
  Full_Batch_Metadata <-
    rbind(Full_Batch_Metadata, Batch1_Meta_SampleID_cleaned)
  Full_Chromatogram_sets <-
    cbind(Full_Chromatogram_sets, Chromatogram_sets)
  
  print(j)
  #this needs to be where the for loop ends
}
#add Scan times to the row names of the full chromatogram set
rownames(Full_Chromatogram_sets) <- rownames(Chromatogram_sets)


#divide intensity values by mass based on SampleBatchID
i <- 1
for (i in 1:length(Full_Batch_Metadata$SampleBatchID)) {
  Sample <- Full_Batch_Metadata$SampleBatchID[i]
  MassVALUE <-
    (Full_Batch_Metadata[which(Full_Batch_Metadata$SampleBatchID == Sample), 8])
  columnofinterest <- which(colnames(Full_Chromatogram_sets) == Sample)
  Full_Chromatogram_sets[, columnofinterest] <-
    Full_Chromatogram_sets[, columnofinterest] / MassVALUE
}


# data is now processed and ready for PCA data is absolute intensity per mg of tissue

write.csv(Full_Batch_Metadata, "Output/Tables/Wild_Helianthus_Metadata.csv")

write.csv(Full_Chromatogram_sets, "Output/Tables/Wild_Helianthus_Chromatograms.csv")



