#RT identification 
library(bclust)



data(gaelle)

gaelle.bclust<-bclust(gaelle,
                      transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
ditplot(gaelle.bclust,varimp=imp(gaelle.bclust)$var,horizbar.plot=TRUE,
        plot.width=5,horizbar.size=0.2,ylab.mar=4)



imp(gaelle.bclust)


View(gaelle)

?imp
#unreplicated clustering

wildtype<-rep(1,55) #initiate a vector
wildtype[c(1:3,48:51,40:43)]<-2 #associate 2 to wildtypes
ditplot(gaelle.bclust,varimp=imp(gaelle.bclust)$var,horizbar.plot=TRUE,
        plot.width=5,horizbar.size=0.2,vertbar=wildtype,
        vertbar.col=c("white","violet"),ylab.mar=4)
#mark wildtype plants using violet 

gaelle.id<-rep(1:14,c(3,rep(4,13))) 
# first 3 rows replication of ColWT, 4 for the rest
gaelle.lab<-c("ColWT","d172","d263","isa2",
              "sex4","dpe2","mex1","sex3","pgm","sex1","WsWT","tpt","RLDWT","ke103")
gaelle.bclust<-bclust(gaelle,rep.id=gaelle.id,
                      labels=gaelle.lab,transformed.par=c(-1.84,-0.99,1.63,0.08,-0.16,-1.68))
ditplot(gaelle.bclust,varimp=imp(gaelle.bclust)$var,horizbar.plot=TRUE)
#replicated clustering







#read in CSV's for metadata
#read in chromatograms transposed so that the sample ID is the ROW and RT is the column
Helianthus_Chromatogram<-t(read.csv("Output/Tables/Wild_Helianthus_Chromatograms.csv",header = T,row.names = 1))
Helianthus_Metadata<-read.csv("Output/Tables/Wild_Helianthus_Metadata.csv",header = T,row.names = 10)
#View(head(Helianthus_Chromatogram))

#remove unreplicated sample
Helianthus_Metadata<-Helianthus_Metadata[!Helianthus_Metadata$Species=="ual(hal)",]

#create column that merges species & treatment
Helianthus_Metadata$COMBO<-factor(paste0(Helianthus_Metadata$Species,"-",Helianthus_Metadata$Treatment))



#merge metatdata and chromatogram for pca

Helianthus_Chromatogram_Metadata<-merge(Helianthus_Metadata,Helianthus_Chromatogram, by=0,all=F)
row.names(Helianthus_Chromatogram_Metadata)<-Helianthus_Chromatogram_Metadata$Row.names
#order based on factor
Helianthus_Chromatogram_Metadata<-Helianthus_Chromatogram_Metadata[order(Helianthus_Chromatogram_Metadata$COMBO),]
dim(Helianthus_Chromatogram_Metadata)
#View(Helianthus_Chromatogram_Metadata[1:20,1:20])

#create dataframes for clusting and hyperparameter estimation
Chromatogram_Rep.ID<-as.integer(Helianthus_Chromatogram_Metadata$COMBO)


Chromatogram.lab<-unique(Helianthus_Chromatogram_Metadata$COMBO)

Chromatogram.bclust<-as.matrix(Helianthus_Chromatogram_Metadata[,12:300])

any(is.na(Chromatogram.bclust))
length(Chromatogram.lab)
View(Chromatogram.bclust[1:20,1:20])
#estimate hyperarameters for use in gaussian model 

meansumsq <- meancss(x=Chromatogram.bclust,rep.id = Chromatogram_Rep.ID)
?meancss
#optimzation function
optimfunc <- function(phi) {-loglikelihood(x.mean = meansumsq$mean, x.css = meansumsq$css, repno = meansumsq$repno, transformed.par = phi)}
x.tpar <- optim(rep(0, 6), optimfunc, method = "BFGS")$par

x.tpar<-loglikelihood(x.mean = meansumsq$mean,x.css = meansumsq$css,repno = meansumsq$repno, transformed.par = rep(0,6))


x.tpar


bclust.obj <- bclust(Chromatogram.bclust, rep.id = Chromatogram_Rep.ID, labels = Chromatogram.lab, transformed.par = x.tpar)



