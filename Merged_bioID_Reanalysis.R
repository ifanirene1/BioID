#Libraries
library(plyr)
library(impute)
library(pcaMethods)
library("preprocessCore")
library("RankProd")
library("sva")
library("RCurl")
setwd("~/Google_Drive/PhD/Project_IF_1/BioID/reanalysis")
#Preprocessing
##merge data
bioid.o9.data <- read.csv("20160620.O9-1.full.csv", header=T)
bioid.3t3.data <- read.csv("3T3_10T_original.csv", header=T)
rownames(bioid.o9.data) <- paste(bioid.o9.data$Gene.ID)
rownames(bioid.3t3.data) <- make.names(bioid.3t3.data$Gene.ID,unique = T)
bioid.o9.data.two.class <-bioid.o9.data[(rowSums(is.na(bioid.o9.data[15:17])) <= 1 | rowSums(is.na(bioid.o9.data[18:20])) <= 1),]
bioid.3t3.data.two.class <-bioid.3t3.data[(rowSums(is.na(bioid.3t3.data[17:19])) <= 1 | rowSums(is.na(bioid.3t3.data[20:22])| rowSums(is.na(bioid.3t3.data[23:24])) <= 1)),]
bioid.merge.data<-merge(bioid.3t3.data.two.class, bioid.o9.data.two.class, by=0, all=TRUE)
rownames(bioid.merge.data) <- paste(bioid.merge.data$Row.names)
bioid.merge.data2<-bioid.merge.data[,c(18:25, 59:64)]
rownames(bioid.merge.data2)<-rownames(bioid.merge.data)



#Sum Normalization
m<-as.matrix(bioid.merge.data2)
bioid.merge.mean<-scale(m, center=FALSE, scale = 5e-05*colSums(m, na.rm = TRUE))
rownames(bioid.merge.mean) <- rownames(bioid.merge.data2)
bioid.merge.log2<-log2(bioid.merge.mean)
ggplot(data=melt(bioid.merge.log2) ,aes(X2,value,colour=X2)) + geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
write.csv(bioid.merge.mean, file = "bioid.merge.mean.csv")
bioid.merge.data.log<-log2(as.matrix(bioid.merge.data2))
ggplot(data=melt(bioid.merge.data.log) ,aes(X2,value,colour=X2)) + geom_boxplot()+theme(axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
boxplot(bioid.merge.data.log, names=c("G1","G2","G3","T1","T2","T3","G1_O9-1","G2_O9-1","G3_O9-1","T1_O9-1","T2_O9-1","T3_O9-1"),ylim= c(-1, 10), las=2)
abline(h=2, col="red", lwd=1)

pc <- pca(t(bioid.merge.log2), center=TRUE, nPcs=12)
pc.scores <- as.data.frame(scores(pc))
#cumulative variance
barplot(pc@R2, ylim=c(0, 0.5))
plot(pc.scores[,1], pc.scores[,2], cex = 2, col = c(rep("green", 3), rep("red", 3),"green","red",rep("green", 3), rep("red", 3)), pch=19,ylim=c(-30, 30),xlim=c(-30, 30))
text(pc.scores[,1], pc.scores[,2], c("G1","G2","G3","T1","T2","T3","G_10T","G_10T",",G1_O9-1","G2_O9-1","G3_O9-1","T1_O9-1","T2_O9-1","T3_O9-1"), cex=0.5, pos=3, col="black")
abline(h = 0, v = 0, col = "gray", lty = 2)
library(dplyr)
setStockcol(NULL)
setStockcol(paste0(getStockcol(), 90))

write.csv(bioid.merge.mean, file="bioid.merge.mean.csv")
row.names(pc@completeObs)
#Analysis

library(MSnbase)
library(msmsTests)
write.csv(bioid.merge.mean, file="bioid.merge.mean.csv")

##Construct feature data

fDataCsv<-read.csv("fData.csv", row.names = 1)
exprsCsv<-read.csv("bioid.merge.mean.csv", row.names = 1)
pDataCsv<-read.csv("pdataFile.csv", row.names = 1)
colnames(exprsCsv)<-row.names(pDataCsv)
res <- readMSnSet(exprsFile="bioid.merge.mean.csv",sep=",", row.names = 1)
#Alternative res <- readMSnSet("bioid.merge.mean.csv", "pdataFile.csv",sep=",", row.names = 1)
sampleNames(res)<-row.names(pDataCsv)
pData(res)<-pDataCsv
#fData(res)<-fDataCsv

###imputation
exp<-pp.msms.data(res)
count<-exprs(exp)
null.f <- "y~1"
alt.f <- "y~treat"
test <- msms.edgeR(exp,alt.f,null.f,div=NULL,fnm="treat")

#Select subset and EdgeR NEGATIVE BINOMIAL
?msms.edgeR
Filter<-pData(exp)$cell=="O9"
o9<- exp[,Filter]
pData(o9)
exp.o9<-exprs(o9)


## for normalisation div <- apply(exprs(exp),2,sum)
null.f <- "y~1"
alt.f <- "y~treat"
test <- msms.edgeR(o9,alt.f,null.f,div=NULL,fnm="treat")
str(test)
head(test)
## values: LogFC-estimated log fold changes, LR-likelihood ratio statistic and corresponding p-value

Filter<-pData(exp)$cell=="3T3"
msn.3t3 <- exp[,Filter]
pData(msn.3t3)
exp.3t3<-exprs(msn.3t3)
View(exp.3t3)
null.f <- "y~1"
alt.f <- "y~treat"


test2 <- msms.edgeR(msn.3t3,alt.f,null.f,div=NULL,fnm="treat")
str(test2)
head(test2)

pData(o9)$treat
nb.tbl.o9 <- test.results(test,o9,pData(o9)$treat,"Twist1","GFP",div=NULL,
                          alpha=0.05,minSpC=2,minLFC=log2(3),method="BH")$tres
nb.tbl.3t3<- test.results(test2,msn.3t3,pData(msn.3t3)$treat,"Twist1","GFP",div=NULL,
                          alpha=0.05,minSpC=2,minLFC=log2(3),method="BH")$tres

?test.results
par(mar=c(5,4,0.5,2)+0.1)
res.volcanoplot(nb.tbl.o9,max.pval=0.05,min.LFC=log2(3),maxx=7,maxy=15,ylbls=5)
res.volcanoplot(nb.tbl.3t3,max.pval=0.05,min.LFC=log2(3),maxx=7,maxy=15,ylbls=3)

?res.volcanoplot()
write.csv(test,file="edgeR.test.O9-1.csv")
write.csv(test2,file="edgeR.test.3T3.csv")
write.csv(nb.tbl, file = "nb.tbl.O9-1.csv")
write.csv(nb.tbl2, file = "nb.tbl.3T3.csv")
write.csv(nb.tbl.o9, file = "nb.tbl.O9-1.2.csv")
write.csv(nb.tbl.3t3, file = "nb.tbl.3T3.2.csv")

pData(res)

####SIMULATION##########

groups<-factor(c(1,1,1,2,2,2,3,4,5,5,5,6,6,6))
m<-count<-exprs(exp)
y<-DGEList(counts= m, group= groups)

groups2<-factor(c(1,1,1,2,2,2,5,5,5,6,6,6))
m2<-m[,c(1:6,9:14)]
y2<-DGEList(counts= m2, group= groups2)
design<-(model.matrix(~0+groups2))
y2<-estimateDisp(y2,design)
fit <- glmQLFit(y2, design)
test.3t3 <- glmQLFTest(fit, contrast=c(-1,1,0,0))
write.table(topTags(test.3t3,n=nrow(y2$counts)), file= "EdgeR_3T3.txt", quote = FALSE, sep="\t",row.names=T)

?glmQLFTest
y2$common.dispersion

test.o9<-glmQLFTest(fit, contrast=c(0,0,-1,1))
write.table(topTags(test.o9,n=nrow(y2$counts)), file= "EdgeR_o9-1.txt", quote = FALSE, sep="\t",row.names=T)

#y= calcNormFactors(y)
plotMDS(b, col=brewer.pal(6,"Set1")[factor(groups)])
design <- model.matrix(~0+groups7, data = y7$samples)
design
y7 <- estimateDisp(y7,design)
fit <- glmQLFit(y7,design)
qlf <- glmQLFTest(fit,coef=7)
topTags<-topTags(qlf)
plotBCV(y7)


m1<-merge[2:19]
rownames(m1)<-merge[,1]
y<-DGEList(counts= m1, group= groups)

y2 <- estimateGLMCommonDisp(y2,design)
y2$common.dispersion^0.5
bcv=0.09981819

bcv <- 0.09981819
#counts <- matrix( rnbinom(40,size=1/bcv^2,mu=10), 20,2)
m3<-m[,7:8]
y3 <- DGEList(counts=m3, group=1:2)
et <- exactTest(y3, dispersion=bcv^2)
et.rank<-topTags(et,n=nrow(y$counts))

write.table(et.rank, file= "10T_dispersion.txt", quote = FALSE, sep="\t",row.names=T)




