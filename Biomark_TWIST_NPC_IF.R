source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")
library(BiocInstaller)
library("HTqPCR")
library(FactoMineR)
library(rafalib)
library(ggrepel)
library(grid)
library(gridExtra)
library(ggplot2)
library(plyr)
library(dplyr)
library(Rtsne)
library(DESeq2)
library(pcaExplorer)

setwd("~/Google_Drive/PhD/Project_IF_1/Biomark")
### FUNCTIONS ############################################################################################################################


#################################
### filtering bad quality data ##
#################################

## The column ``Call'' in the sample file contains information about the result of the qPCR reaction.
## Per default, a call of ``Pass'' is translated into ``OK'' in the  \Rcode{featureCategory},
## and ``Fail'' as ``Undetermined''

filter_data <- function(raw){
  
  ## filter empty and dDris annotated samples
  ## set entire column to NA
  exprs(raw)[,c(grep("H2O", colnames(exprs(raw))),
                grep("-D", colnames(exprs(raw))),
                grep("Water", colnames(exprs(raw))),
                grep("EMPTY", colnames(exprs(raw))),
                grep("Empty", colnames(exprs(raw))),
                grep("-R", colnames(exprs(raw))))]  <- NA
  
  ## ctCall = Failed set to NA
  fail_ind <- which(flag(raw) == "Undetermined")
  exprs(raw)[fail_ind] <- NA
  

  filtered <-  exprs(raw)
  
  ## remove column or row if entire is NA
  filtered_ind_col <- which(colSums(is.na(filtered)) == nrow(filtered))
  if(length(filtered_ind_col) > 0){
    filtered_col <- filtered[,-filtered_ind_col]
  }else{
    filtered_col <- filtered
  }
  
  filtered_ind_row <- which(rowSums(is.na(filtered_col)) == ncol(filtered_col))
  if(length(filtered_ind_row) > 0){
    filtered_row <- filtered_col[-filtered_ind_row,]
  } else{
    filtered_row <- filtered_col
  }
  return(filtered_row)
}



##########################
## Calculate expression ##
##########################

calculate_expression <- function(filtered, LOD){
  exp <- LOD - filtered
  exp[which(exp <= 0)] <- 0.000001
  return(exp)
}



####################
## Remove missing ##
####################

## based on the above assessment have set 40 NA for samples and 23 for genes

remove_missing <- function(exp){
  
  missing <- list()
  
  fail.sample <- which(colSums(is.na(exp))  >= 75)
  if(length(fail.sample) != 0){
    filt.sample <- exp[,-fail.sample]
  }else{
    filt.sample <- exp
  }
  
  ## filter failed genes
  fail.gene <- which(rowSums(is.na(exp))  >= 75)
  if(length(fail.gene) != 0){
    filt.sample2 <- filt.sample[-fail.gene,]
  }else{
    filt.sample2 <- filt.sample
  }
  
  missing <- filt.sample2
  return(missing)
}

####################
## Remove HKGenes ##
####################

remove_HK <- function(exp){
  
  HK <- list()
  
  fail.HK <- which(apply(exp,1,sd, na.rm=T)>= 3)
  if(length(fail.HK) != 0){
    filt.HK <- exp[-fail.HK,]
  }else{
    filt.HK <- exp
  }
  
  HK <- filt.HK
  return(HK)
}

#for future graphics
pd1 = position_dodge(0.2)
pd2 = position_dodge(0.65)

##################
## Read in Data ##############################################################################################################################"
##################

path <-"~/Google_Drive/PhD/Project_IF_1/Biomark"

NPC_TWIST<- "NPC_biomark_20180922.csv"


##open data
NPC_TWIST_samples_path <- paste(path, NPC_TWIST , sep="/")
NPC_TWISTdata <- read.csv("NPC_biomark_20180922.csv", sep=",", header=F)
head(NPC_TWISTdata,n=20)

## remove white spaces
#WNT
NPC_TWISTdata[,2] <- gsub(" ", "", NPC_TWISTdata[,2], fixed = TRUE)


sample_NPC_TWISTdata <- as.vector(unique((NPC_TWISTdata[-c(1:12),2])))

## read in CT data
raw_NPC_TWISTdata <- readCtData(files = NPC_TWIST , format="BioMark", n.features= 96, n.data = 96, samples = sample_NPC_TWISTdata)

dim(raw_NPC_TWISTdata)
rownames(exprs(raw_NPC_TWISTdata))
colnames(exprs(raw_NPC_TWISTdata))


#####################
### filtering data ##
#####################

## doublets triplets, empty etc
filtered_NPC_TWISTdata <- filter_data(raw_NPC_TWISTdata)

###########

### LOD ###
###########

exp_NPC_TWISTdata <- calculate_expression(filtered_NPC_TWISTdata, 28)

####################
## Remove missing ##
####################

NPC_TWISTdata_na_filt <- remove_missing(exp_NPC_TWISTdata)


####################
## Normalisation  ##
####################

#WNT
# Extract the HK genes


NPC_TWISTdata_HK  <- NPC_TWISTdata_na_filt[c(24,36,48),]
NPC_TWISTdata_noHK  <- NPC_TWISTdata_na_filt[-c(24,36,48),]
sd_HK <- apply(NPC_TWISTdata_HK,1,sd, na.rm=T)

NPC_TWISTdata_HK_filt <- remove_HK(NPC_TWISTdata_HK)

# Average the HK genes
NPC_TWISTdata_HK_mean  <- apply (NPC_TWISTdata_HK_filt, 2, mean, na.rm = T)
sd(NPC_TWISTdata_HK_mean)

# Normalisation
NPC_TWISTdata_hk_norm <- sweep(NPC_TWISTdata_noHK, 2, NPC_TWISTdata_HK_mean)
##############
## deltaCt  ##
##############

NPC_TWISTdata_dCT  <- 2^(NPC_TWISTdata_hk_norm)
NPC_TWISTdata_dCT <- 2^(NPC_TWISTdata_na_filt)

#########################
## Select sample ##
#########################

NPC_TWISTdata_dCT <- NPC_TWISTdata_dCT[,-c(grep("NPC", colnames(NPC_TWISTdata_dCT)))]
NPC_TWISTdata_dCT <- NPC_TWISTdata_dCT[,-c(grep("ESC", colnames(NPC_TWISTdata_dCT)))]
NPC_TWISTdata_dCT <- NPC_TWISTdata_dCT[,-c(grep("O9", colnames(NPC_TWISTdata_dCT)))]

#NPC_TWISTdata_dCT <- NPC_TWISTdata_dCT[ ,-which(colnames(NPC_TWISTdata_dCT) %in% c("TC7-1","TC7-2"))]

#########################
## Removed si Samples  ##
#########################



#########################
## Assigned  Samples   ##
#########################

CellLine <- substr(colnames(NPC_TWISTdata_dCT), 1, 6)

CellLine <- as.data.frame(CellLine)
CellLine<- revalue(CellLine$CellLine, c("O9_1_siC_1"="O9_1_siC","O9_1_siC_2"="O9_1_siC","O9_1_siC_3"="O9_1_siC","O9_1_siT_1"="O9_1_siT","O9_1_siT_2"="O9_1_siT","O9_1_siT_3"="O9_1_siT","O9_1_siW_1"="O9_1_siW","O9_1_siW_2"="O9_1_siW","O9_1_siW_3"="O9_1_siW"))
CellLine
CellLine <- as.data.frame(CellLine)
colnames(m)
DIFF <- substr(colnames(m), (nchar(colnames(m))-1), (nchar(colnames(m))-1))
DIFF <- as.data.frame(DIFF)
DIFF <- revalue(DIFF$DIFF, c("E"="ES"))
DIFF <- as.data.frame(DIFF)
DIFF <- revalue(DIFF$DIFF, c("-"="DIFF"))
DIFF <- as.data.frame(DIFF)
DIFF <- revalue(DIFF$DIFF, c("_"="DIFF"))


##########
## PCA  ##
##########

tDat <- t(NPC_TWISTdata_dCT)


# PCA general
mypar(1,2)
PCA_TOT <- PCA(tDat , scale.unit=T,ncp=5, axes = c(1,2))

Type <- CellLine

# Add cell types
DTOT=data.frame(tDat,Type)

# Plot 2d + Barycenter
mypar(2,2)
PCA_TOT = PCA(DTOT, scale.unit=TRUE, ncp=5, quali.sup=ncol(DTOT), graph=T)
x11()
plot.PCA(PCA_TOT, habillage=ncol(DTOT), axes=c(2,3), label=c("none"))

#Extract axis 1 and 2 for individuals (Sample)
PCAcoord <- as.data.frame(PCA_TOT$ind)
PCAcoord12 <- cbind.data.frame(PCAcoord[,1], PCAcoord[,2])

#Extract culture condition
PCAcoord12 <- cbind.data.frame(PCAcoord12, CellLine)

colnames(PCAcoord12) <- c("PC1", "PC2","CellType")


#plot
PCA <- ggplot(PCAcoord12, aes(PC1,PC2, colour= CellType)) +
  geom_point(size=3) +
  #scale_color_manual(values = c("black", "magenta", "darkorange3", "green"))+
  xlab(paste("PC1", "(",round(PCA_TOT$eig[1,2], 2), "% )"))+
  ylab(paste("PC2", "(",round(PCA_TOT$eig[2,2], 2), "% )"))+
  theme_bw()
PCA


##########
## TSNE ##
##########

TSNE_TWIST <- Rtsne(PCAcoord12, perplexity = 1)

TSNEdat <- cbind.data.frame(TSNE_TWIST$Y[,1], TSNE_TWIST$Y[,2], CellLine)

TSNEplot <- ggplot(TSNEdat, aes(TSNE_TWIST$Y[,1], TSNE_TWIST$Y[,2], colour= CellLine))+
  geom_point(aes(shape=CellLine), size = 2)#+xlim(-500,250)+ylim(-200,250)
TSNEplot


#########use pca explorer
write.csv(exprs, "TWIST_NPC_o9.csv")
write.csv(NPC_TWISTdata_dCT, "NPC_hk_dT_2019.csv")
exprs<-exprs(raw_NPC_TWISTdata)
colnames(NPC_TWISTdata_dCT)<-t(CellLine)[,1:37]
as.data.frame(t(CellLine))
head(NPC_TWISTdata_dCT)
class(coldata)
colnames(m7)<-CellLine


#m7<-read.csv(file = "TWISTdata_exp_dT.csv", row.names =1)
m<-NPC_TWISTdata_dCT
head(m)

m[is.na(m)]<-0
count.data<-m#*2000000
count.data<-as.matrix(count.data)
m<-as.data.frame(m)
storage.mode(count.data) = "integer"
coldata <- data.frame(condition=CellLine$CellLine)
rownames(coldata) <- colnames(m)

dds <- DESeqDataSetFromMatrix(countData = count.data, colData = coldata, design = ~ condition)
colData(dds)

dds$condition
dds$condition <- relevel(dds$condition, ref="ESC_A2l")

dds <- DESeq(dds)
pcaExplorer(dds=dds)

estimateSizeFactors(dds)
sizeFactors(dds) <- rep(c(1), each=24)
res <- results(dds)
colData(dds)$sizeFactor


head(assay(dds))

head(counts(dds,normalized=TRUE))
########################subsetting##########################
m<-NPC_TWISTdata_dCT
m <- subset(m , select=-c(NPC_AC7_3,EB_AC7_2,grep("O9", colnames(m))))
#m <- subset(m , select=-c(NPC_AC7_3,EB_AC7_2,grep("NPC", colnames(m)))
#drops <- c("x","z")
#DF[ , !(names(DF) %in% drops)]
CellLine <- substr(colnames(m), 1, 7)
CellLine <- as.data.frame(CellLine)
CellLine

Subset<-dds$condition!="NPC_T8-"
dds<-dds[,Subset]
dds$condition<-droplevels(dds$condition)
dds <- estimateSizeFactors(dds)
idx <- rowSums( counts(dds, normalized=TRUE) >= 5 ) >= 3
dds <- dds[idx,]
dds <- DESeq(dds)





m1<-as.matrix(m)
m2<-m1/m1[,72]
m1<-scale(t(m1))
