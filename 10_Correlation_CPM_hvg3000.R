
#correlation analysis
library(corrplot)
library(psych)
library(ggplot2)
Ourstudy <- read.table("Ourstudy_pigCPM3000.txt",sep="\t",header=T,row.names = 1)

dim(Ourstudy)
mg_OurWang <- read.table("mg_OurWang_pigCPM3000.txt",sep="\t",header=T,row.names = 1)

dim(mg_OurWang)

identity_name  <- intersect(row.names(Ourstudy),row.names(mg_OurWang))
com_Ourstudy   <- Ourstudy[identity_name,]
com_mg_OurWang <- mg_OurWang[identity_name,]

#----------------------------------------------------------------------------------------
#pearson
cor_pearson <- cor(as.matrix(com_Ourstudy),as.matrix(com_mg_OurWang),method = "pearson")
write.table(cor_pearson,file="corr_pearson.txt",quote = F,sep="\t")

pdf("corrplot_pearson.pdf",width = 15,height=13)
corrplot(cor_pearson,
         method="circle",
         type="full",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100),
         bg="white",
         is.corr=F,
         add=F,
         diag=T,
         tl.srt=90, 
         cl.pos="r", 
         tl.col="black")
dev.off()      

#----------------------------------------------------------------------------------------
#spearman
cor_spearman <- cor(as.matrix(com_Ourstudy),as.matrix(com_mg_OurWang),method = "spearman")
write.table(aa,file="corr_spearman.txt",quote = F,sep="\t")

pdf("corrplot_spearman.pdf",width = 18,height=16)
corrplot(cor_spearman,
         method="circle",
         type="full",
         col=colorRampPalette(c('#0000ff','#ffffff','#ff0000'))(100),
         bg="white",
         is.corr=F,
         add=F,
         diag=T,
         tl.srt=90, 
         cl.pos="r", 
         tl.col="black")
dev.off()     


#----------------------------------------------------------------------------------------
#pseudo correlation

library(ggplot2)
library(dplyr)
library(ggpubr)
com_Ourstudy   <- Ourstudy[identity_name,]
com_mg_OurWang <- mg_OurWang[identity_name,]
rowMean_com_Ourstudy  <-  rowMeans(com_Ourstudy)
rowMean_com_mg_OurWang <-  rowMeans(com_mg_OurWang)
data=as.data.frame(cbind(rowMean_com_Ourstudy,rowMean_com_mg_OurWang))
colnames(data)=c("rowMean_com_Ourstudy_value","rowMean_com_mg_OurWang_value")

cor(rowMean_com_Ourstudy,rowMean_com_mg_OurWang,method="pearson")
pdf("corrplot_Mean_pearson.pdf",width = 8,height=6)
ggplot(data,aes(x=rowMean_com_Ourstudy,y=rowMean_com_mg_OurWang))+ 
      theme_bw()+
      geom_point(size=1,shape=15)+
      geom_smooth(method=lm)+
      stat_cor(method = 'pearson')
dev.off()


pdf("corrplot_Mean_spearman.pdf",width = 8,height=6)
ggplot(data,aes(x=rowMean_com_Ourstudy,y=rowMean_com_mg_OurWang))+ 
       theme_bw()+
       geom_point(size=1,shape=15)+
       geom_smooth(method=lm)+
       stat_cor(method = 'spearman')
dev.off()
























