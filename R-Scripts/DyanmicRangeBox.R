library("dynRB")
library(dynRB)
library(ggplot2)
library(reshape2)
library(vegan)
library(RColorBrewer)
library(magrittr)
library(BayesFactor)
library(EnhancedVolcano)
library(gridExtra)
library(sm)
library(factoextra)
source("R-Scripts/JordanVolcano.R")
#read in data to make a more useful metadata file and rename 






#chemical data

Data<- read.csv("Data/Voc.csv")

#View(Data)
#standardize data by sample weight

 Data2<-t(sweep(as.matrix(t(Data[c(6:ncol(Data))])), 2, as.vector(Data$weight), `/`))

 #keep for total peak area
 Data[c(6:ncol(Data))]<-Data2

 
 #convert to relative ratio or ignore to use ion counts
#Data2<-t(Data2)
#Data3<-t(prop.table(Data2,margin = 2))


#replace data vaules
#Data[c(6:ncol(Data))]<-Data3
#Data[c(6:ncol(Data))]<-Data2[,-1]

Data$Species<-as.factor(Data$Species)
Data$treatment<-as.factor(Data$treatment)


#Global PCA


Chromatogram.PCA<- prcomp(Data[,6:ncol(Data)],center = T, scale. = T)

#assess eigen structure





#assess individual loadings


PCA.IND <-
  fviz_pca_ind(
    Chromatogram.PCA,
    axes = c(1, 2),
    geom.ind = "point",
    pointshape = 19,
    #label = Data$Species,
    col.ind = Data$treatment,
    mean.point = F
  )


#PCA.IND
#row.names(PCA.IND$data)<-PCA.IND$data$Col.
PCADATA<-merge(PCA.IND$data,Data$Species, by=0, all=T)
#View(PCADATA)
pdf("Output/Figures/Helianthus_plots.pdf")
fviz_eig(Chromatogram.PCA, choice = "variance", addlabels = T)

ggplot(
  PCADATA,
  aes(
    x = x,
    y = y.x,
    col = factor(y.y)
    
  )
) + geom_point(aes(shape=factor(Col.))) + xlab(PCA.IND$labels$x)+ylab(PCA.IND$labels$y)+
  labs(col= "Species",shape= "Treatment")

# plot densities
sm.density.compare(PCADATA$y.x,PCADATA$Col., xlab="PC1")
title(main=bquote(italic("H. annuus")))
# add legend
colfill<-c(2:(2+length(levels(PCADATA$Col.))))
legend('topright',levels(PCADATA$Col.), fill=colfill)

#create species specific analyses
dev.off()



#generate volcano plots 
#subset data according to treatment
CONTROL<-Data[Data$treatment=="Control",]
VOLATILE<-Data[Data$treatment=="Volatile",]
EATEN<-Data[Data$treatment=="Herbivory",]

#empty dataframe to store species specific results
BayesfactorEvidence<-data.frame()
for (j in 6:ncol(Data)) {
  
  
  #calculate volatile difference
  if(sum(VOLATILE[,j])!=0 && sum(CONTROL[,j])!=0){
    #Bayesian factor analysis 
    volatile_bf<-ttestBF(x=VOLATILE[,j],y=CONTROL[,j], paired = FALSE)
    #pull out evidence it needs to be saved as a vector to be paired with LogFC
    volatile_Evidence<-as.data.frame(volatile_bf)[,1:2]
    rownames(volatile_Evidence)<-c()
    #calculate Log2FC
    #replace 0s with 0.00001
    x1cmpd<-VOLATILE[,j]
    x1cmpd[x1cmpd==0]<-0.000001
    x2cmpd<-CONTROL[,j]
    x2cmpd[x2cmpd==0]<-0.000001
    FC<- mean(x1cmpd)/mean(x2cmpd)
    VOLATILElog2FC<-log(FC,2)
    volatile_Eff<-((mean(x1cmpd)-mean(x2cmpd))/sd(c(x1cmpd,x2cmpd)))
    
    
    #combine evidence and log2FC
    Volatile_Vectortoadd<-as.data.frame(cbind(volatile_Evidence,VOLATILElog2FC,volatile_Eff))
    colnames(Volatile_Vectortoadd)<-c("volatile_BF","Volatile_bferror","VOLATILElog2FC","volatile_Eff")
    
  }else{
    Volatile_Vectortoadd<-as.data.frame(cbind(0,0,0,0))
    colnames(Volatile_Vectortoadd)<-c("volatile_BF","Volatile_bferror","VOLATILElog2FC","volatile_Eff")
  }
  
  
  #calculate eaten difference
  if(sum(EATEN[,j])!=0 && sum(CONTROL[,j])!=0){
    #Bayesian factor analysis 
    EATEN_bf<-ttestBF(x=EATEN[,j],y=CONTROL[,j], paired = FALSE)
    #pull out evidence it needs to be saved as a vector to be paired with LogFC
    EATEN_Evidence<-as.data.frame(EATEN_bf)[,1:2]
    rownames(EATEN_Evidence)<-c()
    #calculate Log2FC
    #replace 0s with 0.00001
    x1cmpd<-EATEN[,j]
    x1cmpd[x1cmpd==0]<-0.000001
    x2cmpd<-CONTROL[,j]
    x2cmpd[x2cmpd==0]<-0.000001
    FC<- mean(x1cmpd)/mean(x2cmpd)
    EATENlog2FC<-log(FC,2)
    EATEN_Eff<-((mean(x1cmpd)-mean(x2cmpd))/sd(c(x1cmpd,x2cmpd)))
    
    
    
    #combine evidence and log2FC
    EATEN_Vectortoadd<-as.data.frame(cbind(EATEN_Evidence,EATENlog2FC,EATEN_Eff))
    colnames(EATEN_Vectortoadd)<-c("EATEN_BF","EATEN_bferror","EATENlog2FC","EATEN_Eff")
  }else{
    EATEN_Vectortoadd<-as.data.frame(cbind(0,0,0,0))
    colnames(EATEN_Vectortoadd)<-c("EATEN_BF","EATEN_bferror","EATENlog2FC","EATEN_Eff")
    
  }
  
  
  #calculate volatile_eaten difference
  if(sum(EATEN[,j])!=0 && sum(VOLATILE[,j])!=0){
    #Bayesian factor analysis 
    volatile_eaten_bf<-ttestBF(x=EATEN[,j],y=VOLATILE[,j], paired = FALSE)
    #pull out evidence it needs to be saved as a vector to be paired with LogFC
    volatile_eaten_Evidence<-as.data.frame(volatile_eaten_bf)[,1:2]
    rownames(volatile_eaten_Evidence)<-c()
    #calculate Log2FC
    #replace 0s with 0.00001
    x1cmpd<-EATEN[,j]
    x1cmpd[x1cmpd==0]<-0.000001
    x2cmpd<-VOLATILE[,j]
    x2cmpd[x2cmpd==0]<-0.000001
    FC<- mean(x1cmpd)/mean(x2cmpd)
    volatile_eatenlog2FC<-log(FC,2)
    
    FC<- mean(x1cmpd)/mean(x2cmpd)
    volatile_eatenlog2FC<-log(FC,2)
    volatile_eaten_Eff<-((mean(x1cmpd)-mean(x2cmpd))/sd(c(x1cmpd,x2cmpd)))
    
    
    #combine evidence and log2FC
    volatile_eaten_Vectortoadd<-as.data.frame(cbind(volatile_eaten_Evidence,volatile_eatenlog2FC,volatile_eaten_Eff))
    colnames(volatile_eaten_Vectortoadd)<-c("volatile_eaten_BF","volatile_eaten_bferror","volatile_eatenlog2FC","volatile_eaten_Eff")
  }else{
    volatile_eaten_Vectortoadd<-as.data.frame(cbind(0,0,0,0))
    colnames(volatile_eaten_Vectortoadd)<-c("volatile_eaten_BF","volatile_eaten_bferror","volatile_eatenlog2FC","volatile_eaten_Eff")
  }
  
  #combine data sets
  
  
  
  BayesfactorEvidence<-rbind(BayesfactorEvidence,cbind(Volatile_Vectortoadd,EATEN_Vectortoadd,volatile_eaten_Vectortoadd))
  
}



#end of chemistry loop
#add chemistry names to rows
row.names(BayesfactorEvidence)<-colnames(Data[,6:ncol(Data)])

#write dataframe to directory

write.csv(BayesfactorEvidence,paste0("Ouput/Figures/","/Helianthus_BayesfactorEvidence.csv"))
#here is where we would want to start making volcanoplots

#replace BF 0s with 0.00001

BayesfactorEvidence[BayesfactorEvidence==0]<-0.000001



FCvolatileplot<-JordanVolcano(BayesfactorEvidence,
                              lab = rownames(BayesfactorEvidence),
                              title =bquote(italic("H. annuus") ~ 'volatile-induced vs  control') ,
                              subtitle = "Peak area",
                              FCcutoff = 1,
                              pCutoff = 1,
                              boxedLabels = F,
                              drawConnectors = F,
                              ylim = c(0.001,(max(BayesfactorEvidence$volatile_BF)+1)),
                              xlim = c(-(max(abs(BayesfactorEvidence$VOLATILElog2FC))+1),(max(abs(BayesfactorEvidence$VOLATILElog2FC))+1)),
                              x = 'VOLATILElog2FC',
                              ylab = "Bayes Factor",
                              y = 'volatile_BF')
EFFvolatileplot<-JordanVolcano(BayesfactorEvidence,
                               lab = rownames(BayesfactorEvidence),
                               title =bquote(italic("H. annuus") ~ 'volatile-induced vs  control') ,
                               subtitle = "Peak area",
                               FCcutoff = 0.2,
                               pCutoff = 1,
                               boxedLabels = F,
                               drawConnectors = F,
                               ylim = c(0.001,(max(BayesfactorEvidence$volatile_BF)+1)),
                               xlim = c(-(max(ceiling(abs(BayesfactorEvidence$volatile_Eff)))),max(ceiling(abs(BayesfactorEvidence$volatile_Eff)))),
                               xlab = paste0("Cohen's D"),
                               x = 'volatile_Eff',
                               ylab = "Bayes Factor",
                               legendLabels = c('NS', expression(italic(D)),
                                                'BF', expression(BF~and~italic(D))),
                               y = 'volatile_BF')


eatenplot<-JordanVolcano(BayesfactorEvidence,
                         lab = rownames(BayesfactorEvidence),
                         title =bquote(italic("H. annuus") ~ 'herbivory-induced vs control') ,
                         subtitle = "Peak area",
                         FCcutoff = 0.5,
                         pCutoff = 1,
                         boxedLabels = F,
                         drawConnectors = F,
                         ylim = c(0.001,(max(BayesfactorEvidence$EATEN_BF)+1)),
                         xlim = c(-(max(abs(BayesfactorEvidence$EATENlog2FC))+1),(max(abs(BayesfactorEvidence$EATENlog2FC))+1)),
                         x = 'EATENlog2FC',
                         ylab = "Bayes Factor",
                         y = 'EATEN_BF')
EFFeatenplot<-JordanVolcano(BayesfactorEvidence,
                            lab = rownames(BayesfactorEvidence),
                            title =bquote(italic("H. annuus") ~ 'herbivory-induced vs  control') ,
                            subtitle = "Peak area",
                            FCcutoff = 0.2,
                            pCutoff = 1,
                            boxedLabels = F,
                            drawConnectors = F,
                            ylim = c(0.001,(max(BayesfactorEvidence$EATEN_BF)+1)),
                            xlim = c(-(max(ceiling(abs(BayesfactorEvidence$EATEN_Eff)))),max(ceiling(abs(BayesfactorEvidence$EATEN_Eff)))),
                            xlab = paste0("Cohen's D"),
                            x = 'EATEN_Eff',
                            ylab = "Bayes Factor",
                            legendLabels = c('NS', expression(italic(D)),
                                             'BF', expression(BF~and~italic(D))),
                            y = 'EATEN_BF')


volatileeatenplot<-JordanVolcano(BayesfactorEvidence,
                                 lab = rownames(BayesfactorEvidence),
                                 title =bquote(italic("H. annuus") ~ 'herbivory vs volatile') ,
                                 subtitle = "Peak area",
                                 FCcutoff = 1,
                                 pCutoff = 1,
                                 boxedLabels = F,
                                 drawConnectors = F,
                                 ylim = c(0.001,(max(BayesfactorEvidence$volatile_eaten_BF)+1)),
                                 xlim = c(-(max(abs(BayesfactorEvidence$volatile_eatenlog2FC))+1),(max(abs(BayesfactorEvidence$volatile_eatenlog2FC))+1)),
                                 x = 'volatile_eatenlog2FC',
                                 ylab = "Bayes Factor",
                                 y = 'volatile_eaten_BF')

EFFvolatileeatenplot<-JordanVolcano(BayesfactorEvidence,
                                    lab = rownames(BayesfactorEvidence),
                                    title =bquote(italic("H. annuus") ~ 'herbivory vs  control') ,
                                    subtitle = "Peak area",
                                    FCcutoff = 0.2,
                                    pCutoff = 1,
                                    boxedLabels = F,
                                    drawConnectors = F,
                                    ylim = c(0.001,(max(BayesfactorEvidence$volatile_eaten_BF)+1)),
                                    xlim = c(-(max(ceiling(abs(BayesfactorEvidence$volatile_eaten_Eff)))),max(ceiling(abs(BayesfactorEvidence$volatile_eaten_Eff)))),
                                    xlab = paste0("Cohen's D"),
                                    x = 'volatile_eaten_Eff',
                                    ylab = "Bayes Factor",
                                    legendLabels = c('NS', expression(italic(D)),
                                                     'BF', expression(BF~and~italic(D))),
                                    y = 'volatile_eaten_BF')





pdf(paste0("Ouput/Figures/","HELIANTHUS_VOLCANO_plots.pdf"))
print(FCvolatileplot)
print(EFFvolatileplot)
print(eatenplot)
print(EFFeatenplot)
print(volatileeatenplot)
print(EFFvolatileeatenplot)
dev.off()











#get list of species names
SpeciesNames<-unique(Data$Species)

i<-1



for (i in 1:length(SpeciesNames)) {
  print(paste(SpeciesNames[i],"Starting..."))
#create new dataframe based on evidence and Logfold change

  
  
  
  
  

#create folder to store data 
dir.create(path=paste0("Ouput/Figures/",SpeciesNames[i],"/"),recursive = T)

#species data subset
SpeciesData<-Data[Data$Species==SpeciesNames[i],]


#PCA at species level
#create dataset with only compounds with observations
# 
# PCA_SpeciesData<-SpeciesData[,6:ncol(SpeciesData)]
# 
# PCA_SpeciesData<- PCA_SpeciesData[,-(which(colSums(PCA_SpeciesData)==0))]
# SpeciesData[,6:ncol(SpeciesData)]<-PCA_SpeciesData
# #PCA
# Chromatogram.PCA<- prcomp(SpeciesData[,6:ncol(SpeciesData)],center = T, scale. = T)
# 
# #assess eigen structure
# 
# 
# fviz_eig(Chromatogram.PCA, choice = "variance", addlabels = T)
# 
# 
# #assess individual loadings
# PCA.IND <-
#   fviz_pca_ind(
#     Chromatogram.PCA,
#     axes = c(1, 2),
#     geom.ind = "point",
#     pointshape = 19,
#     #label = Data$Species,
#     col.ind = SpeciesData$treatment,
#     mean.point = F
#   )
# 
# 
# PCA.IND
# 
# #PCA.IND
# PCADATA<-PCA.IND$data
# 
# # plot densities
# sm.density.compare(PCADATA$x,PCADATA$Col., xlab="PC1")
# title(main=bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i])))))
# # add legend
# colfill<-c(2:(2+length(levels(PCADATA$Col.))))
# legend('topright',levels(PCADATA$Col.), fill=colfill)




#subset data according to treatment
CONTROL<-SpeciesData[SpeciesData$treatment=="Control",]
VOLATILE<-SpeciesData[SpeciesData$treatment=="Volatile",]
EATEN<-SpeciesData[SpeciesData$treatment=="Herbivory",]



#empty dataframe to store species specific results
BayesfactorEvidence<-data.frame()
for (j in 6:ncol(Data)) {
  

#calculate volatile difference
if(sum(VOLATILE[,j])!=0 && sum(CONTROL[,j])!=0){
  #Bayesian factor analysis 
  volatile_bf<-ttestBF(x=VOLATILE[,j],y=CONTROL[,j], paired = FALSE)
  #pull out evidence it needs to be saved as a vector to be paired with LogFC
  volatile_Evidence<-as.data.frame(volatile_bf)[,1:2]
  rownames(volatile_Evidence)<-c()
  #calculate Log2FC
  #replace 0s with 0.00001
  x1cmpd<-VOLATILE[,j]
  #x1cmpd[x1cmpd==0]<-0.000000000000000000001
  x2cmpd<-CONTROL[,j]
  #x2cmpd[x2cmpd==0]<-0.000000000000000000001
  FC<- mean(x1cmpd)/mean(x2cmpd)
  VOLATILElog2FC<-log(FC,2)
  volatile_Eff<-((mean(x1cmpd)-mean(x2cmpd))/sd(c(x1cmpd,x2cmpd)))
  
  
  
  
  #combine evidence and log2FC
  Volatile_Vectortoadd<-data.frame()
  Volatile_Vectortoadd<-as.data.frame(cbind(volatile_Evidence,VOLATILElog2FC,volatile_Eff))
  colnames(Volatile_Vectortoadd)<-c("volatile_BF","Volatile_bferror","VOLATILElog2FC","volatile_Eff")
  
}else{
  Volatile_Vectortoadd<-data.frame()
  Volatile_Vectortoadd<-as.data.frame(cbind(0,0,0,0))
  colnames(Volatile_Vectortoadd)<-c("volatile_BF","Volatile_bferror","VOLATILElog2FC","volatile_Eff")
}


#calculate eaten difference
if(sum(EATEN[,j])!=0 && sum(CONTROL[,j])!=0){
  #Bayesian factor analysis 
  EATEN_bf<-ttestBF(x=EATEN[,j],y=CONTROL[,j], paired = FALSE)
  #pull out evidence it needs to be saved as a vector to be paired with LogFC
  EATEN_Evidence<-as.data.frame(EATEN_bf)[,1:2]
  rownames(EATEN_Evidence)<-c()
  #calculate Log2FC
  #replace 0s with 0.00001
  x1cmpd<-EATEN[,j]
  #x1cmpd[x1cmpd==0]<-0.000000000000000000001
  x2cmpd<-CONTROL[,j]
  #x2cmpd[x2cmpd==0]<-0.000000000000000000001
  FC<- mean(x1cmpd)/mean(x2cmpd)
  EATENlog2FC<-log(FC,2)
  EATEN_Eff<-((mean(x1cmpd)-mean(x2cmpd))/sd(c(x1cmpd,x2cmpd)))
  
  
  
  #combine evidence and log2FC
  EATEN_Vectortoadd<-data.frame()
  EATEN_Vectortoadd<-as.data.frame(cbind(EATEN_Evidence,EATENlog2FC,EATEN_Eff))
  colnames(EATEN_Vectortoadd)<-c("EATEN_BF","EATEN_bferror","EATENlog2FC","EATEN_Eff")
}else{
  EATEN_Vectortoadd<-data.frame()
  EATEN_Vectortoadd<-as.data.frame(cbind(0,0,0,0))
  colnames(EATEN_Vectortoadd)<-c("EATEN_BF","EATEN_bferror","EATENlog2FC","EATEN_Eff")
  
  }


#calculate volatile_eaten difference
if(sum(EATEN[,j])!=0 && sum(VOLATILE[,j])!=0){
  #Bayesian factor analysis 
  volatile_eaten_bf<-ttestBF(x=EATEN[,j],y=VOLATILE[,j], paired = FALSE)
  #pull out evidence it needs to be saved as a vector to be paired with LogFC
  volatile_eaten_Evidence<-as.data.frame(volatile_eaten_bf)[,1:2]
  rownames(volatile_eaten_Evidence)<-c()
  #calculate Log2FC
  #replace 0s with 0.00001
  x1cmpd<-EATEN[,j]
  #x1cmpd[x1cmpd==0]<-0.000000000000000000001
  x2cmpd<-VOLATILE[,j]
  #x2cmpd[x2cmpd==0]<-0.000000000000000000001
  FC<- mean(x1cmpd)/mean(x2cmpd)
  volatile_eatenlog2FC<-log(FC,2)
  
  FC<- mean(x1cmpd)/mean(x2cmpd)
  volatile_eatenlog2FC<-log(FC,2)
  volatile_eaten_Eff<-((mean(x1cmpd)-mean(x2cmpd))/sd(c(x1cmpd,x2cmpd)))
  
  
  #combine evidence and log2FC
  volatile_eaten_Vectortoadd<-data.frame()
  volatile_eaten_Vectortoadd<-as.data.frame(cbind(volatile_eaten_Evidence,volatile_eatenlog2FC,volatile_eaten_Eff))
  colnames(volatile_eaten_Vectortoadd)<-c("volatile_eaten_BF","volatile_eaten_bferror","volatile_eatenlog2FC","volatile_eaten_Eff")
}else{
  volatile_eaten_Vectortoadd<-data.frame()
  volatile_eaten_Vectortoadd<-as.data.frame(cbind(0,0,0,0))
colnames(volatile_eaten_Vectortoadd)<-c("volatile_eaten_BF","volatile_eaten_bferror","volatile_eatenlog2FC","volatile_eaten_Eff")
  }

#combine data sets



BayesfactorEvidence<-rbind(BayesfactorEvidence,cbind(Volatile_Vectortoadd,EATEN_Vectortoadd,volatile_eaten_Vectortoadd))

}



#end of chemistry loop
#add chemistry names to rows
row.names(BayesfactorEvidence)<-colnames(Data[,6:ncol(Data)])

#write dataframe to directory

write.csv(BayesfactorEvidence,paste0("Ouput/Figures/",SpeciesNames[i],"/",SpeciesNames[i],"_BayesfactorEvidence.csv"))
#here is where we would want to start making volcanoplots

#replace BF 0s with 0.00001

BayesfactorEvidence[BayesfactorEvidence==0]<-0.000000000000000000001

FCvolatileplot<-JordanVolcano(BayesfactorEvidence,
                              lab = rownames(BayesfactorEvidence),
                            title =bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i]))) ~ 'volatile-induced vs  control') ,
                            subtitle = "Peak area",
                              FCcutoff = 1,
                              pCutoff = 1,
                              boxedLabels = F,
                              drawConnectors = F,
                              ylim = c(0.001,(max(BayesfactorEvidence$volatile_BF)+1)),
                              xlim = c(-(max(abs(BayesfactorEvidence$VOLATILElog2FC))+1),(max(abs(BayesfactorEvidence$VOLATILElog2FC))+1)),
                              x = 'VOLATILElog2FC',
                              ylab = "Bayes Factor",
                              y = 'volatile_BF')
EFFvolatileplot<-JordanVolcano(BayesfactorEvidence,
                              lab = rownames(BayesfactorEvidence),
                              title =bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i]))) ~ 'volatile-induced vs  control') ,
                              subtitle = "Peak area",
                              FCcutoff = 0.2,
                              pCutoff = 1,
                              boxedLabels = F,
                              drawConnectors = F,
                              ylim = c(0.001,(max(BayesfactorEvidence$volatile_BF)+1)),
                              xlim = c(-(max(ceiling(abs(BayesfactorEvidence$volatile_Eff)))),max(ceiling(abs(BayesfactorEvidence$volatile_Eff)))),
                              xlab = paste0("Cohen's D"),
                              x = 'volatile_Eff',
                              ylab = "Bayes Factor",
                              legendLabels = c('NS', expression(italic(D)),
                                               'BF', expression(BF~and~italic(D))),
                              y = 'volatile_BF')


eatenplot<-JordanVolcano(BayesfactorEvidence,
                           lab = rownames(BayesfactorEvidence),
                         title =bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i]))) ~ 'herbivory-induced vs control') ,
                         subtitle = "Peak area",
                         FCcutoff = 0.5,
                           pCutoff = 1,
                           boxedLabels = F,
                           drawConnectors = F,
                           ylim = c(0.001,(max(BayesfactorEvidence$EATEN_BF)+1)),
                           xlim = c(-(max(abs(BayesfactorEvidence$EATENlog2FC))+1),(max(abs(BayesfactorEvidence$EATENlog2FC))+1)),
                           x = 'EATENlog2FC',
                           ylab = "Bayes Factor",
                           y = 'EATEN_BF')
EFFeatenplot<-JordanVolcano(BayesfactorEvidence,
                               lab = rownames(BayesfactorEvidence),
                               title =bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i]))) ~ 'herbivory-induced vs  control') ,
                               subtitle = "Peak area",
                               FCcutoff = 0.2,
                               pCutoff = 1,
                               boxedLabels = F,
                               drawConnectors = F,
                               ylim = c(0.001,(max(BayesfactorEvidence$EATEN_BF)+1)),
                               xlim = c(-(max(ceiling(abs(BayesfactorEvidence$EATEN_Eff)))),max(ceiling(abs(BayesfactorEvidence$EATEN_Eff)))),
                               xlab = paste0("Cohen's D"),
                               x = 'EATEN_Eff',
                               ylab = "Bayes Factor",
                               legendLabels = c('NS', expression(italic(D)),
                                                'BF', expression(BF~and~italic(D))),
                               y = 'EATEN_BF')


volatileeatenplot<-JordanVolcano(BayesfactorEvidence,
                                   lab = rownames(BayesfactorEvidence),
                                 title =bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i]))) ~ 'herbivory vs volatile') ,
                                 subtitle = "Peak area",
                                 FCcutoff = 1,
                                   pCutoff = 1,
                                   boxedLabels = F,
                                   drawConnectors = F,
                                   ylim = c(0.001,(max(BayesfactorEvidence$volatile_eaten_BF)+1)),
                                   xlim = c(-(max(abs(BayesfactorEvidence$volatile_eatenlog2FC))+1),(max(abs(BayesfactorEvidence$volatile_eatenlog2FC))+1)),
                                   x = 'volatile_eatenlog2FC',
                                   ylab = "Bayes Factor",
                                   y = 'volatile_eaten_BF')

EFFvolatileeatenplot<-JordanVolcano(BayesfactorEvidence,
                            lab = rownames(BayesfactorEvidence),
                            title =bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i]))) ~ 'herbivory vs volatile') ,
                            subtitle = "Peak area",
                            FCcutoff = 0.2,
                            pCutoff = 1,
                            boxedLabels = F,
                            drawConnectors = F,
                            ylim = c(0.001,(max(BayesfactorEvidence$volatile_eaten_BF)+1)),
                            xlim = c(-(max(ceiling(abs(BayesfactorEvidence$volatile_eaten_Eff)))),max(ceiling(abs(BayesfactorEvidence$volatile_eaten_Eff)))),
                            xlab = paste0("Cohen's D"),
                            x = 'volatile_eaten_Eff',
                            ylab = "Bayes Factor",
                            legendLabels = c('NS', expression(italic(D)),
                                             'BF', expression(BF~and~italic(D))),
                            y = 'volatile_eaten_BF')





pdf(paste0("Ouput/Figures/",SpeciesNames[i],"/",SpeciesNames[i],"_plots.pdf"))
#print(PCA.IND)
# plot densities
#sm.density.compare(PCADATA$x,PCADATA$Col., xlab="PC1")
#title(main=bquote(italic("H.") ~ italic(.(paste(SpeciesNames[i])))))
# add legend
#colfill<-c(2:(2+length(levels(PCADATA$Col.))))
#legend('topright',levels(PCADATA$Col.), fill=colfill)
print(ggarrange(FCvolatileplot, EFFvolatileplot, eatenplot, EFFeatenplot, 
                ncol = 2, nrow = 2, labels = c("A.","B.","C.","D.","E.","F."), align = "hv"))

print(FCvolatileplot)
print(EFFvolatileplot)
print(eatenplot)
print(EFFeatenplot)
print(volatileeatenplot)
print(EFFvolatileeatenplot)


dev.off()





#end Species for loop here
print(paste(SpeciesNames[i],"Finished..."))
}






#add new stuff here









#save chains to species folder with plots
chains <- posterior(bf, iterations = 100000)
summary(chains)
plot(chains)




# Estimate preprocessing parameters
preproc.param <- Data3 %>% 
  preProcess(method = c("center", "scale"))
# Transform the data using the estimated parameters
train.transformed <- preproc.param %>% predict(Data3)

heatmap(as.matrix(train.transformed[1:24,]))
model <- lda(Data$treatment~Data3)
View(model$means)
wine.lda.values <- predict(model)

wine.lda.values$posterior

plot(wine.lda.values$x[,1],wine.lda.values$x[,2],col=as.factor(Data$Species)) # make a scatterplot
text(wine.lda.values$x[,1],wine.lda.values$x[,2],Data$Species,cex=0.7,pos=4,col="red") # add labels


plot.lda(model)

View(model$)


View(scatterplotMatrix(Data3[1:2,]))


#volcano plot

Groupmeans<-t(model$means)


eFC<- Groupmeans[,2]/Groupmeans[,3]

View(eFC)

log2FC<-log(eFC,2)

pvalue<-apply(log2FC,1,function(x){t.test(Groupmeans[,2],Groupmeans[,1])})



View(log2FC)



Chromatogram.PCA<-prcomp(train.transformed[1:12,])

fviz_eig(Chromatogram.PCA, choice = "variance", addlabels = T)

PCA.IND <-
  fviz_pca_ind(
    Chromatogram.PCA,
    axes = c(3, 2),
    geom.ind = "point",
    pointshape = 19,
    #label = SAMMETA$CORE12,
    col.ind = Data$treatment[1:12],
    mean.point = F
  )
PCA.IND

corrplot::corrplot(cor(Data3))

unique(Data$Species)



#Dynamic range box examination

#convert 0 to missing data only works with missing data if theres still one value 
#so skip or adjust

#Ready.Data[Ready.Data == 0] <- NA

#add metadata to first column

dynb.ready.data<-data.frame(Data$treatment,Data3)
dynb.ready.data$Data.treatment %<>% factor




#overall overlap
HypervolumeOverlap <- dynRB_VPa(dynb.ready.data)

#overlap by dimension
DimensionOverlap <- dynRB_Pn(dynb.ready.data)

#quantify dimension size
DimensionSize <- dynRB_Vn(dynb.ready.data)

View(HypervolumeOverlap$result)





