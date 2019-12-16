if(!require(devtools)) install.packages("devtools")#install the package for significance annotation, good for publication quality figures
devtools::install_github("kassambara/ggpubr")
library(plyr)
library(dplyr)
library(ggplot2)
library(ggpubr)
setwd("~/Google_Drive/PhD/Project_IF_1/Biomark")
Dat<-read.csv(file = "IWP2_TW_dCT.csv", row.names =1)
Match <- c("ES_A", "T2.3.", "T8.2.","A2loxcre")#select cell types for comparison
Dat2<-Dat[,c(grep(c(paste(Match,collapse="|")),colnames(Dat)))]
head(Dat2)
colnames(Dat2)
#Match <- c("siC_", "siT_", "siC8","siTC8")#select cell types for comparison
#Dat2<-Dat[,c(grep(c(paste(Match,collapse="|")),colnames(Dat)))]

CellType <- substr(colnames(Dat2), 1, 3)

CellType <- as.data.frame(CellType)
CellType<- revalue(CellType$CellType, c("ES_"="ESC","A2l"="WT","T2."="Twist1+/-","T8."="Twist1-/-"))#fix or rename them
CellType <- as.data.frame(CellType)
CellType

#colnames(Dat2)<-t(CellType)


myplots1 <- list()


#This loop will generate one plot per gene and save a pdf file with each plot on a separate page
for (i in 1:length(rownames(Dat2))) {
  
  n <- i
  #normalisation
  Dat3 <- Dat2[n,]
  norm <- min(Dat3[1,], na.rm=T)
  Dat4 <- t(Dat3/norm)
 Dat4[is.na(Dat4)]<-0.000000001
  Dat4<-cbind.data.frame(Dat4,CellType)
  colnames(Dat4)<-c("Value", "Cell")
  #comparisons <- list(c("siControl","siT"), c("siControl","siChd8"),c("siControl","siTC8"))#set groups for stat comparison, or use ref.group = "WT" as argument in stat_compare_means()
  p<-ggplot(Dat4, aes(Cell,Value))+
    stat_summary(fun.y=mean,colour="red", geom="point")+
    geom_dotplot (aes(fill=Cell, colour=Cell), binaxis = "y", dotsize= 1, stackdir = "center")+
    stat_summary(fun.data = mean_se, geom = "errorbar", colour= "red", width= 0.5)+
    ggtitle(rownames(Dat[n,]))+
    theme_bw()+
    theme(aspect.ratio = 1, text = element_text(size=20),axis.text.x = element_text(angle=90, hjust=1))+
    stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group = "WT")+ # Add pairwise comparisons p-value
    stat_compare_means(method= "anova")+theme(legend.position="none")     # Add global p-value,default is wilcox
  myplots1[[i]] <- p
}

pdf("twist1_null_EB_mannWhitney.pdf")
for (i in 1:length(rownames(Dat2))) {
  print(myplots1[[i]])
}
dev.off()


D<-t(Dat2)
D[is.na(D)]<-0.000000001
head(D)
D<-as.matrix(D)
D1<-t(D[3,])
D<-sweep(D,2,D1,FUN = "/")

Dat5<-as.data.frame(D)
Dat5<-cbind.data.frame(Celltype=CellType,Dat5)
head(Dat5)
Dat5$CellType<-factor(Dat5$CellType, levels =c("ESC","WT","Twist1+/-","Twist1-/-") )
ggboxplot(Dat5, x="CellType",y = c("Sox2","Sox9","Pou5f1"), 
          merge = "flip",
          color = "Celltype",width = 0.7,
          ylab = "Log2FC", add = "jitter", add.params = list(size = 0.1, jitter = 0.2))#+position_dodge(0.6)#+
  #stat_compare_means(label = "p.signif", method = "t.test",comparisons = list(c("ESC", "EB_"),c("ESC", "NPC")))
?stat_compare_means
?ggboxplot
ggbox



anno_df = compare_means(c("Sox2","Sox9","Pou5f1") ~ CellType, group.by = "CellType", data = Dat5) %>%
  mutate(y_pos = 40)


ggplot(Dat2, aes(x=Celltype, y = c("Sox2","Sox9","Pou5F1"))) + 
  geom_boxplot(position=position_dodge()) + 
  geom_point(aes(color=supp), position=position_jitterdodge()) + 
  facet_wrap(~Celltype) + 
  ggsignif::geom_signif(
    data=anno_df, 
    aes(xmin=group1, xmax=group2, annotations=p.adj, y_position=y_pos), 
    manual=TRUE
  )



ggbarplot(Dat5, x="CellType",y = c("Sox2","Sox9","Pou5f1"), 
          merge = "flip",
          color = "CellType",width = 0.7,
          ylab = "Log2FC", add = "jitter", add.params = list(size = 0.1, jitter = 0.2))+
stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group ="WT", paired=FALSE)#comparisons = list(c("ESC", "EB_"),c("ESC", "NPC")))


ggbarplot(Dat5, x="CellType",y = c("Wnt3","Fzd10","Axin2","Prickle1"), 
          merge = "flip",
          width = 0.7,add = "mean_se",
          ylab = "Fold Change")
 
  stat_compare_means(aes(group = CellType),label = "p.signif", method = "wilcox")#, comparisons = list(c("WT", "Twist1+/-"),c("WT", "Twist1-/-")))


  c("Wnt3","Fzd10","Axin2")
  c("Gsc","Pdx1","Rbm47","Foxa2","Gata6","Sox17")
  c("Fgf10","Zic1","Zeb2","Sox9","FoxD4")
  c("Dppa2","Sox2","Nanog","Klf4","Pou3f1","Pou5f1")
  "Flk1","Mesp1","T/Bra","Mixl1","Pdgfra","Vegfa","Snail2","Cdh1"
  
  
  
  
  
  
  ##########################Normality test####################
  na.omit(Dat5$Foxa2)
  
  shapiro.test(log2(na.omit(Dat5$Flk1)))
  
  shapiro.test(Dat6$Snail2)
  
  shapiro.test(Dat6$Klf4)
  shapiro.test(Dat5[,12])
Dat6<-Dat5[grep("ES",Dat5$CellType),]





ggbarplot(Dat5, x="CellType",y = "Sox2", 
          merge = "flip",
          color = "CellType",width = 0.7,
          ylab = "Log2FC", add = "jitter", add.params = list(size = 0.1, jitter = 0.2))+
  stat_compare_means(label = "p.signif", method = "wilcox.test",ref.group ="WT")#comparisons = list(c("ESC", "EB_"),c("ESC", "NPC")))



