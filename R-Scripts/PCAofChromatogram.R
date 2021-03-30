#PCA of chromatogram data


#packages used

library(factoextra)
library(missMDA)
library(corrplot)
library(ggfortify)
library(BayesFactor)
library(ggpubr)
library(plyr)

#read in CSV's for metadata
#read in chromatograms transposed so that the sample ID is the ROW and RT is the column
Helianthus_Chromatogram<-t(read.csv("Output/Tables/Wild_Helianthus_Chromatograms.csv",header = T,row.names = 1))
Helianthus_Metadata<-read.csv("Output/Tables/Wild_Helianthus_Metadata.csv",header = T,row.names = 10)
View(head(Helianthus_Metadata))



#summary table SpeciesXTreatment
count_freq <- count(Helianthus_Metadata, c("Species", "Treatment"))
View(count_freq)
#merge metatdata and chromatogram for pca

Helianthus_Chromatogram_Metadata<-merge(Helianthus_Metadata,Helianthus_Chromatogram, by=0,all=T)
row.names(Helianthus_Chromatogram_Metadata)<-Helianthus_Chromatogram_Metadata$Row.names
#PCA

Chromatogram.PCA<- prcomp(Helianthus_Chromatogram_Metadata[,11:ncol(Helianthus_Chromatogram_Metadata)],center = T, scale. = T)


#assess eigen structure


fviz_eig(Chromatogram.PCA, choice = "variance", addlabels = T)

View(Chromatogram.PCA$rotation)
#assess individual loadings

PCA.IND <-
  fviz_pca_ind(
    Chromatogram.PCA,
    axes = c(1, 2),
    geom.ind = "point",
    pointshape = 19,
    #label = SAMMETA$CORE12,
    col.ind = Helianthus_Chromatogram_Metadata$Row.names,
    mean.point = F
  )

#PCA.IND
row.names(PCA.IND$data)<-PCA.IND$data$Col.
PCADATA<-merge(PCA.IND$data,Helianthus_Chromatogram_Metadata[,7:8], by=0, all=T)


ggplot(
  PCADATA,
  aes(
    x = x,
    y = y,
    col = factor(Species)
  )
) + geom_point(aes(shape=factor(Treatment)))# + scale_shape_manual(values = c(21,22,24))
  #scale_size_manual(values = c(5, 2)) +
  
  #ggtitle(PCA.IND$labels$title) + 
  #xlab(PCA.IND$labels$x) + ylab(PCA.IND$labels$y) +
  #labs(fill = "Breeding group", shape = "Core 12") + 
  #guides(size = FALSE) + geom_hline(yintercept =0,linetype = "dashed",color = "black") + geom_vline(xintercept = 0,linetype = "dashed",color = "black") + theme_bw()


PCADATA_Angustifolous<-PCADATA[PCADATA$Species=="agrestis",]
#View(PCADATA_Angustifolous)

View(PCADATA)

ggplot(
  PCADATA_Angustifolous,
  aes(
    x = x,
    y = y,
    col = factor(Treatment)
  )
) + geom_point(aes(shape=factor(Treatment)))






#species level PCA



Helianthus_Chromatogram_Metadata_Species<-Helianthus_Chromatogram_Metadata[Helianthus_Chromatogram_Metadata$Species=="silphioides",]

#View(Helianthus_Chromatogram_Metadata_Species)

#PCA

Chromatogram.PCA<- prcomp(Helianthus_Chromatogram_Metadata_Species[,11:ncol(Helianthus_Chromatogram_Metadata)],center = T, scale. = T)

View(Chromatogram.PCA$rotation)
#assess eigen structure


fviz_eig(Chromatogram.PCA, choice = "variance", addlabels = T)


#assess individual loadings

PCA.IND <-
  fviz_pca_ind(
    Chromatogram.PCA,
    axes = c(5, 6),
    geom.ind = "point",
    pointshape = 19,
    #label = SAMMETA$CORE12,
    col.ind = Helianthus_Chromatogram_Metadata_Species$Row.names,
    mean.point = F
  )

#PCA.IND
row.names(PCA.IND$data)<-PCA.IND$data$Col.
PCADATA<-merge(PCA.IND$data,Helianthus_Chromatogram_Metadata_Species[,7:8], by=0, all=T)

ggplot(
  PCADATA,
  aes(
    x = x,
    y = y,
    col = factor(Treatment)
  )
) + geom_point(aes(shape=factor(Treatment)))


View(PCADATA)
