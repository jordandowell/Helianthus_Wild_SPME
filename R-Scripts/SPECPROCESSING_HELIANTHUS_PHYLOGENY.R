####SPECTRAL PREDICTION OF CONTEXT
#packages used
library(readtext)
library(dplyr)
library(tidyr)
library(pls)
library(devtools)
devtools::install_github("richardjtelford/ggbiplot", ref = "experimental")
library(ggbiplot)
library(tidyverse)
library(ggpubr)
library(ggplot2) # visualization
library(gridExtra)
library(e1071)
library(caret)
#read all files in folder

# set workind directory-->setwd("~/filefolderpath")



file_list <- list.files(path = "Data/HyperspecReflectance/",pattern = "*.txt",full.names = T)
#extract metadata
metadata<- readtext::readtext("Data/HyperspecReflectance/*.txt",
                              docvarsfrom = "filenames", 
                              docvarnames = c("species", "treatment", "plant_replicate", "measurement", "measurement_replicate"),
                              dvsep = "_", 
                              encoding = NULL)

length(unique(metadata$species))
#extract measurements
#setwd("~/Documents/Documents - jordan’s MacBook Pro (2)/Scans/Ground Tissue (Raw Data)/ALLscandsgroundADD")
# update this file path to point toward appropriate folder on your computer
#blank matrix
dataset <- as.data.frame(matrix(nrow = 1, ncol = 2554))
for (file in file_list) {
  # if the file name contains reflectance spliced remove 14row header
  if (grepl("reflectancespliced", file) == T) {
    v1 <- (read.table(
      file,
      header = F,
      sep = "\t",
      skip = 14
    ))
    v1 <- t(v1[-1])
    dataset <- rbind(dataset, v1)
  }
  # if the file name does not containts reflectance spliced remove 14row header
  if (grepl("reflectancespliced", file) == F) {
    v2 <- (read.table(
      file,
      header = F,
      sep = "\t",
      skip = 14
    ))
    v2 <- t(v2[-1])
    dataset <- rbind(dataset, v2)
  }
}

View(dataset[1:10,1:10])
dataset <- dataset[-1, ]
#extract wavelengths
rows <- (read.table(
  file_list[1],
  header = F,
  sep = "\t",
  skip = 14
))
rows <- t(rows[-2])
colnames(dataset) <- rows
#bind data
fulldata <- cbind(metadata, dataset)
dim(fulldata)
View(fulldata[1:10,1:10])



write.csv(fulldata, "allscansNOTAVERAGED_07.12.19.csv")
#i saved the big dataset
#LArgerdataset<-fulldata

#for loop for Species PCA figures 

for (j in 1:length(unique(fulldata$species))) {
  print(unique(fulldata$species)[j])
#subsetdata to species 
  
  SpeciesData_PCA<- fulldata[fulldata$species==unique(fulldata$species)[j],]

  
  # # create PCA to reduce colinearity
  # #method is singular value decomposition
  SpeciesData_PCA.pca <- prcomp(SpeciesData_PCA[, 8:length(SpeciesData_PCA)],
                       center = T,
                       scale. = T)
  # construct plots
 
  onetwo<-autoplot(SpeciesData_PCA.pca,x=1,y=2, data = SpeciesData_PCA, colour ='treatment')+ 
    geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
    geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
    theme_bw()
  onethree<-autoplot(SpeciesData_PCA.pca,x=1,y=3, data = SpeciesData_PCA, colour ='treatment')+ 
    geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
    geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
    theme_bw()
  onefour<-autoplot(SpeciesData_PCA.pca,x=1,y=4, data = SpeciesData_PCA, colour ='treatment')+ 
    geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
    geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
    theme_bw()
  twothree<-autoplot(SpeciesData_PCA.pca,x=2,y=3, data = SpeciesData_PCA, colour ='treatment')+ 
    geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
    geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
    theme_bw()
  twofour<-autoplot(SpeciesData_PCA.pca,x=2,y=4, data = SpeciesData_PCA, colour ='treatment')+ 
    geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
    geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
    theme_bw()
  threefour<-autoplot(SpeciesData_PCA.pca,x=3,y=4, data = SpeciesData_PCA, colour ='treatment')+ 
    geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
    geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
    theme_bw()
  
  
  pdf(paste0("Output/",unique(fulldata$species)[j],"_PCAplots.pdf"),onefile = FALSE)
  print(ggarrange(onetwo, onethree, onefour, twothree,twofour,threefour, 
            ncol = 2, nrow = 3, labels = c("A.","B.","C.","D.","E.","F."), align = "hv",common.legend = T))
  dev.off()

  }







#All classes
#filter Test & training Data for Global model
#find testing and training data 
Metacombos<- cbind(metadata$species,metadata$treatment,metadata$plant_replicate)
Metacombos %>% group_by_all %>% count
res <- Metacombos[!duplicated(Metacombos), ]
ind<- duplicated(res[,1:2])
testingreplicates<-as.data.frame(res[!ind,])

colnames(testingreplicates)<-c("species","treatment","plant_replicate")
testingreplicates[,3]<-as.integer(testingreplicates[,3])

#get training replicates

res <- Metacombos[!duplicated(Metacombos), ]
ind<- duplicated(res[,1:2])
trainingreplicates<-as.data.frame(res[ind,])
colnames(trainingreplicates)<-c("species","treatment","plant_replicate")
trainingreplicates[,3]<-as.integer(trainingreplicates[,3])


#67% of data
Global.Training <- match_df(fulldata,trainingreplicates)
#33% of data
Global.Test <- match_df(fulldata,testingreplicates)
View(Global.Test[1:5,1:30])
View(fulldata[1:5,1:30])
#averaged <- as.data.frame(colMeans(fulldata[,8:ncol(fulldata)]))
#write.csv(averaged, "maxaverage.csv")

dim(fulldata)
#only control and treatment
#CEfulldata <- fulldata %>% filter(treatement %in% c("e","c"))


# # create PCA to reduce colinearity
# #method is singular value decomposition
 Global.pca <- prcomp(Global.Training[, 8:length(Global.Training)],
                      center = T,
                      scale. = T)
# #assess >= 95% variance explanation
summary(Global.pca)

autoplot(Global.pca,x=1,y=4, data = Global.Training, colour ='species', shape="treatment")+ 
  geom_hline(yintercept =0,linetype = "dashed",color = "black") + 
  geom_vline(xintercept = 0,linetype = "dashed",color = "black") + 
  theme_bw()




fviz_eig(Global.pca, choice = "variance", addlabels = T)
#PCA 1 = 52.11%; 2 = 23.2%; 3 = 11.42%; 4 = 8.49%
#project Training data into PCA space
Global.Training.Projection <- predict(Global.pca, Global.Training)
#bind first 8 components to new Training dataset
Global.Training.PCA <-
  data.frame(Global.Training[, 1:7], Global.Training.Projection[, 1:4])
#project Test data into PCA space
Global.Test.Projection <- predict(Global.pca, Global.Test)
#bind first 4 components to new Test dataset
Global.Test.PCA <-
 data.frame(Global.Test[, 1:7], Global.Test.Projection[, 1:4])









#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
Global.svm <-
  tune.svm(
    y = factor(Global.Training.PCA$treatment),
    x = Global.Training.PCA[,8:ncol(Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    #data = Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.Global.SVM <-
  svm(
    y = factor(Global.Training.PCA$treatment),
    x = Global.Training.PCA[,8:ncol(Global.Training.PCA)],
    kernel = "radial",
    cost = Global.svm$best.parameters$cost,
    gamma = Global.svm$best.parameters$gamma,
    #data = Global.Training.PCA, 
    cross = 10,
    probability = TRUE)
  
length(unique(Global.Training.PCA$species))
summary(Tuned.Global.SVM)
#Select tuned model and create confusion matrix of Test data
Global.Training.Prediction <-
  predict(Tuned.Global.SVM,
          Global.Training.PCA[, 8:ncol(Global.Training.PCA)],
          probability = TRUE)
Global.Training.Confusion <-
  confusionMatrix(factor(Global.Training.Prediction),
                  factor(Global.Training.PCA$treatment))
Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
Global.Test.Prediction <-
  predict(Tuned.Global.SVM,
          Global.Test.PCA[, 8:ncol(Global.Test.PCA)],
          probability = TRUE)
Global.Test.Confusion <-
  confusionMatrix(factor(Global.Test.Prediction),
                  factor(Global.Test.PCA$treatment))
Global.Test.Confusion



##Senarios 1-6
SENARIO1<- paste(c("maxa", "gig","gros","nut"), collapse = "|")
SENARIO2<- paste(c("atr", "silph", "ang", "mol"), collapse = "|")
SENARIO3<- paste(c("agro", "pet","debi", "gra"), collapse = "|")
SENARIO4<- paste(c("agro", "maxa", "atr", "ariz"), collapse = "|")
SENARIO5<- paste(c("pet", "gig", "ang", "lac"), collapse = "|")
SENARIO6<- paste(c("debi", "gros", "silph", "mol"), collapse = "|")

#########SENARIO 1

#filter Test & training Data for Global model
# #67% of data
S1.Global.Training <- fulldata %>% filter(!grepl(SENARIO1, species))
# #33% of data
S1.Global.Test <- fulldata %>% filter(grepl(SENARIO1, species))


# 
# # create PCA from training data to reduce collinearity
# #method is singular value decomposition
S1.Global.pca <- prcomp(S1.Global.Training[, 8:length(S1.Global.Training)],
                     center = T,
                     scale. = T)
# #assess >= 95% variance explanation
summary(S1.Global.pca)
#PCA 1 = 51%; 2 = 24%; 3 = 12%; 4 = 8%
#project Training data into PCA space
S1.Global.Training.Projection <- predict(S1.Global.pca, S1.Global.Training)
#bind first 4 components to new Training dataset
S1.Global.Training.PCA <-
  data.frame(S1.Global.Training[, 1:7], S1.Global.Training.Projection[, 1:4])
#project Test data into PCA space
S1.Global.Test.Projection <- predict(S1.Global.pca, S1.Global.Test)
#bind first 4 components to new Test dataset
S1.Global.Test.PCA <-
  data.frame(S1.Global.Test[, 1:7], S1.Global.Test.Projection[, 1:4])



#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
S1.Global.svm <-
  tune.svm(
    y = factor(S1.Global.Training.PCA$treatement),
    x = S1.Global.Training.PCA[,8:ncol(S1.Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    data = S1.Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.S1.Global.SVM <-
  svm(
    y = factor(S1.Global.Training.PCA$treatement),
    x = S1.Global.Training.PCA[,8:ncol(S1.Global.Training.PCA)],
    kernel = "radial",
    cost = S1.Global.svm$best.parameters$cost,
    gamma = S1.Global.svm$best.parameters$gamma,
    data = S1.Global.Training.PCA, 
    cross = 10,
    probability = TRUE)


summary(Tuned.S1.Global.SVM)
#Select tuned model and create confusion matrix of Test data
S1.Global.Training.Prediction <-
  predict(Tuned.S1.Global.SVM,
          S1.Global.Training.PCA[, 8:ncol(S1.Global.Training.PCA)],
          probability = TRUE)
S1.Global.Training.Confusion <-
  confusionMatrix(factor(S1.Global.Training.Prediction),
                  factor(S1.Global.Training.PCA$treatement))
S1.Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
S1.Global.Test.Prediction <-
  predict(Tuned.S1.Global.SVM,
          S1.Global.Test.PCA[, 8:ncol(S1.Global.Test.PCA)],
          probability = TRUE)
S1.Global.Test.Confusion <-
  confusionMatrix(factor(S1.Global.Test.Prediction),
                  factor(S1.Global.Test.PCA$treatement))
S1.Global.Test.Confusion



#########SENARIO 2

#filter Test & training Data for Global model
# #67% of data
S2.Global.Training <- fulldata %>% filter(!grepl(SENARIO2, species))
# #33% of data
S2.Global.Test <- fulldata %>% filter(grepl(SENARIO2, species))


# 
# # create PCA from training data to reduce collinearity
# #method is singular value decomposition
S2.Global.pca <- prcomp(S2.Global.Training[, 8:length(S2.Global.Training)],
                        center = T,
                        scale. = T)
# #assess >= 95% variance explanation
summary(S2.Global.pca)
#PCA 1 = 51%; 2 = 24%; 3 = 12%; 4 = 8%
#project Training data into PCA space
S2.Global.Training.Projection <- predict(S2.Global.pca, S2.Global.Training)
#bind first 4 components to new Training dataset
S2.Global.Training.PCA <-
  data.frame(S2.Global.Training[, 1:7], S2.Global.Training.Projection[, 1:4])
#project Test data into PCA space
S2.Global.Test.Projection <- predict(S2.Global.pca, S2.Global.Test)
#bind first 4 components to new Test dataset
S2.Global.Test.PCA <-
  data.frame(S2.Global.Test[, 1:7], S2.Global.Test.Projection[, 1:4])



#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
S2.Global.svm <-
  tune.svm(
    y = factor(S2.Global.Training.PCA$treatement),
    x = S2.Global.Training.PCA[,8:ncol(S2.Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    data = S2.Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.S2.Global.SVM <-
  svm(
    y = factor(S2.Global.Training.PCA$treatement),
    x = S2.Global.Training.PCA[,8:ncol(S2.Global.Training.PCA)],
    kernel = "radial",
    cost = S2.Global.svm$best.parameters$cost,
    gamma = S2.Global.svm$best.parameters$gamma,
    data = S2.Global.Training.PCA, 
    cross = 10,
    probability = TRUE)


summary(Tuned.S2.Global.SVM)
#Select tuned model and create confusion matrix of Test data
S2.Global.Training.Prediction <-
  predict(Tuned.S2.Global.SVM,
          S2.Global.Training.PCA[, 8:ncol(S2.Global.Training.PCA)],
          probability = TRUE)
S2.Global.Training.Confusion <-
  confusionMatrix(factor(S2.Global.Training.Prediction),
                  factor(S2.Global.Training.PCA$treatement))
S2.Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
S2.Global.Test.Prediction <-
  predict(Tuned.S2.Global.SVM,
          S2.Global.Test.PCA[, 8:ncol(S2.Global.Test.PCA)],
          probability = TRUE)
S2.Global.Test.Confusion <-
  confusionMatrix(factor(S2.Global.Test.Prediction),
                  factor(S2.Global.Test.PCA$treatement))
S2.Global.Test.Confusion


#########SENARIO 3

#filter Test & training Data for Global model
# #67% of data
S3.Global.Training <- fulldata %>% filter(!grepl(SENARIO3, species))
# #33% of data
S3.Global.Test <- fulldata %>% filter(grepl(SENARIO3, species))


# 
# # create PCA from training data to reduce collinearity
# #method is singular value decomposition
S3.Global.pca <- prcomp(S3.Global.Training[, 8:length(S3.Global.Training)],
                        center = T,
                        scale. = T)
# #assess >= 95% variance explanation
summary(S3.Global.pca)
#PCA 1 = 51%; 2 = 24%; 3 = 12%; 4 = 8%
#project Training data into PCA space
S3.Global.Training.Projection <- predict(S3.Global.pca, S3.Global.Training)
#bind first 4 components to new Training dataset
S3.Global.Training.PCA <-
  data.frame(S3.Global.Training[, 1:7], S3.Global.Training.Projection[, 1:4])
#project Test data into PCA space
S3.Global.Test.Projection <- predict(S3.Global.pca, S3.Global.Test)
#bind first 4 components to new Test dataset
S3.Global.Test.PCA <-
  data.frame(S3.Global.Test[, 1:7], S3.Global.Test.Projection[, 1:4])



#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
S3.Global.svm <-
  tune.svm(
    y = factor(S3.Global.Training.PCA$treatement),
    x = S3.Global.Training.PCA[,8:ncol(S3.Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    data = S3.Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.S3.Global.SVM <-
  svm(
    y = factor(S3.Global.Training.PCA$treatement),
    x = S3.Global.Training.PCA[,8:ncol(S3.Global.Training.PCA)],
    kernel = "radial",
    cost = S3.Global.svm$best.parameters$cost,
    gamma = S3.Global.svm$best.parameters$gamma,
    data = S3.Global.Training.PCA, 
    cross = 10,
    probability = TRUE)


summary(Tuned.S3.Global.SVM)
#Select tuned model and create confusion matrix of Test data
S3.Global.Training.Prediction <-
  predict(Tuned.S3.Global.SVM,
          S3.Global.Training.PCA[, 8:ncol(S3.Global.Training.PCA)],
          probability = TRUE)
S3.Global.Training.Confusion <-
  confusionMatrix(factor(S3.Global.Training.Prediction),
                  factor(S3.Global.Training.PCA$treatement))
S3.Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
S3.Global.Test.Prediction <-
  predict(Tuned.S3.Global.SVM,
          S3.Global.Test.PCA[, 8:ncol(S3.Global.Test.PCA)],
          probability = TRUE)
S3.Global.Test.Confusion <-
  confusionMatrix(factor(S3.Global.Test.Prediction),
                  factor(S3.Global.Test.PCA$treatement))
S3.Global.Test.Confusion



#########SENARIO 4

#filter Test & training Data for Global model
# #67% of data
S4.Global.Training <- fulldata %>% filter(!grepl(SENARIO4, species))
# #33% of data
S4.Global.Test <- fulldata %>% filter(grepl(SENARIO4, species))


# 
# # create PCA from training data to reduce collinearity
# #method is singular value decomposition
S4.Global.pca <- prcomp(S4.Global.Training[, 8:length(S4.Global.Training)],
                        center = T,
                        scale. = T)
# #assess >= 95% variance explanation
summary(S4.Global.pca)
#PCA 1 = 51%; 2 = 24%; 3 = 12%; 4 = 8%
#project Training data into PCA space
S4.Global.Training.Projection <- predict(S4.Global.pca, S4.Global.Training)
#bind first 4 components to new Training dataset
S4.Global.Training.PCA <-
  data.frame(S4.Global.Training[, 1:7], S4.Global.Training.Projection[, 1:4])
#project Test data into PCA space
S4.Global.Test.Projection <- predict(S4.Global.pca, S4.Global.Test)
#bind first 4 components to new Test dataset
S4.Global.Test.PCA <-
  data.frame(S4.Global.Test[, 1:7], S4.Global.Test.Projection[, 1:4])



#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
S4.Global.svm <-
  tune.svm(
    y = factor(S4.Global.Training.PCA$treatement),
    x = S4.Global.Training.PCA[,8:ncol(S4.Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    data = S4.Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.S4.Global.SVM <-
  svm(
    y = factor(S4.Global.Training.PCA$treatement),
    x = S4.Global.Training.PCA[,8:ncol(S4.Global.Training.PCA)],
    kernel = "radial",
    cost = S4.Global.svm$best.parameters$cost,
    gamma = S4.Global.svm$best.parameters$gamma,
    data = S4.Global.Training.PCA, 
    cross = 10,
    probability = TRUE)


summary(Tuned.S4.Global.SVM)
#Select tuned model and create confusion matrix of Test data
S4.Global.Training.Prediction <-
  predict(Tuned.S4.Global.SVM,
          S4.Global.Training.PCA[, 8:ncol(S4.Global.Training.PCA)],
          probability = TRUE)
S4.Global.Training.Confusion <-
  confusionMatrix(factor(S4.Global.Training.Prediction),
                  factor(S4.Global.Training.PCA$treatement))
S4.Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
S4.Global.Test.Prediction <-
  predict(Tuned.S4.Global.SVM,
          S4.Global.Test.PCA[, 8:ncol(S4.Global.Test.PCA)],
          probability = TRUE)
S4.Global.Test.Confusion <-
  confusionMatrix(factor(S4.Global.Test.Prediction),
                  factor(S4.Global.Test.PCA$treatement))
S4.Global.Test.Confusion




#########SENARIO 5

#filter Test & training Data for Global model
# #67% of data
S5.Global.Training <- fulldata %>% filter(!grepl(SENARIO5, species))
# #33% of data
S5.Global.Test <- fulldata %>% filter(grepl(SENARIO5, species))


# 
# # create PCA from training data to reduce collinearity
# #method is singular value decomposition
S5.Global.pca <- prcomp(S5.Global.Training[, 8:length(S5.Global.Training)],
                        center = T,
                        scale. = T)
# #assess >= 95% variance explanation
summary(S5.Global.pca)
#PCA 1 = 51%; 2 = 24%; 3 = 12%; 4 = 8%
#project Training data into PCA space
S5.Global.Training.Projection <- predict(S5.Global.pca, S5.Global.Training)
#bind first 4 components to new Training dataset
S5.Global.Training.PCA <-
  data.frame(S5.Global.Training[, 1:7], S5.Global.Training.Projection[, 1:4])
#project Test data into PCA space
S5.Global.Test.Projection <- predict(S5.Global.pca, S5.Global.Test)
#bind first 4 components to new Test dataset
S5.Global.Test.PCA <-
  data.frame(S5.Global.Test[, 1:7], S5.Global.Test.Projection[, 1:4])



#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
S5.Global.svm <-
  tune.svm(
    y = factor(S5.Global.Training.PCA$treatement),
    x = S5.Global.Training.PCA[,8:ncol(S5.Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    data = S5.Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.S5.Global.SVM <-
  svm(
    y = factor(S5.Global.Training.PCA$treatement),
    x = S5.Global.Training.PCA[,8:ncol(S5.Global.Training.PCA)],
    kernel = "radial",
    cost = S5.Global.svm$best.parameters$cost,
    gamma = S5.Global.svm$best.parameters$gamma,
    data = S5.Global.Training.PCA, 
    cross = 10,
    probability = TRUE)


summary(Tuned.S5.Global.SVM)
#Select tuned model and create confusion matrix of Test data
S5.Global.Training.Prediction <-
  predict(Tuned.S5.Global.SVM,
          S5.Global.Training.PCA[, 8:ncol(S5.Global.Training.PCA)],
          probability = TRUE)
S5.Global.Training.Confusion <-
  confusionMatrix(factor(S5.Global.Training.Prediction),
                  factor(S5.Global.Training.PCA$treatement))
S5.Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
S5.Global.Test.Prediction <-
  predict(Tuned.S5.Global.SVM,
          S5.Global.Test.PCA[, 8:ncol(S5.Global.Test.PCA)],
          probability = TRUE)
S5.Global.Test.Confusion <-
  confusionMatrix(factor(S5.Global.Test.Prediction),
                  factor(S5.Global.Test.PCA$treatement))
S5.Global.Test.Confusion



#########SENARIO 6

#filter Test & training Data for Global model
# #67% of data
S6.Global.Training <- fulldata %>% filter(!grepl(SENARIO6, species))
# #33% of data
S6.Global.Test <- fulldata %>% filter(grepl(SENARIO6, species))


# 
# # create PCA from training data to reduce collinearity
# #method is singular value decomposition
S6.Global.pca <- prcomp(S6.Global.Training[, 8:length(S6.Global.Training)],
                        center = T,
                        scale. = T)
# #assess >= 95% variance explanation
summary(S6.Global.pca)
#PCA 1 = 51%; 2 = 24%; 3 = 12%; 4 = 8%
#project Training data into PCA space
S6.Global.Training.Projection <- predict(S6.Global.pca, S6.Global.Training)
#bind first 4 components to new Training dataset
S6.Global.Training.PCA <-
  data.frame(S6.Global.Training[, 1:7], S6.Global.Training.Projection[, 1:4])
#project Test data into PCA space
S6.Global.Test.Projection <- predict(S6.Global.pca, S6.Global.Test)
#bind first 4 components to new Test dataset
S6.Global.Test.PCA <-
  data.frame(S6.Global.Test[, 1:7], S6.Global.Test.Projection[, 1:4])



#Begin SVM modeling


#parameter tuning
#logarithmic parameter search
#LOOCV cross valdiation 
#Radial basis kernel 
S6.Global.svm <-
  tune.svm(
    y = factor(S6.Global.Training.PCA$treatement),
    x = S6.Global.Training.PCA[,8:ncol(S6.Global.Training.PCA)],
    kernel = "radial",
    cost = 10 ^ (-5:4),
    gamma = 10 ^ (-5:4),
    data = S6.Global.Training.PCA,
    tunecontrol = tune.control(cross = 10, nrepeat = 10)
  )

#create tuned model with optimal parameters
#repeated LOOCV ~approximate same model
#set probabilities = to true to get class likelihood
Tuned.S6.Global.SVM <-
  svm(
    y = factor(S6.Global.Training.PCA$treatement),
    x = S6.Global.Training.PCA[,8:ncol(S6.Global.Training.PCA)],
    kernel = "radial",
    cost = S6.Global.svm$best.parameters$cost,
    gamma = S6.Global.svm$best.parameters$gamma,
    data = S6.Global.Training.PCA, 
    cross = 10,
    probability = TRUE)


summary(Tuned.S6.Global.SVM)
#Select tuned model and create confusion matrix of Test data
S6.Global.Training.Prediction <-
  predict(Tuned.S6.Global.SVM,
          S6.Global.Training.PCA[, 8:ncol(S6.Global.Training.PCA)],
          probability = TRUE)
S6.Global.Training.Confusion <-
  confusionMatrix(factor(S6.Global.Training.Prediction),
                  factor(S6.Global.Training.PCA$treatement))
S6.Global.Training.Confusion
#1-sensitivity = False control
#1-specificty = False induction
#1-accuracy is OOB

#Select tuned model and create confusion matrix of Test data
S6.Global.Test.Prediction <-
  predict(Tuned.S6.Global.SVM,
          S6.Global.Test.PCA[, 8:ncol(S6.Global.Test.PCA)],
          probability = TRUE)
S6.Global.Test.Confusion <-
  confusionMatrix(factor(S6.Global.Test.Prediction),
                  factor(S6.Global.Test.PCA$treatement))
S6.Global.Test.Confusion
#FINISHED
