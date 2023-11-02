library(clusterProfiler)
library(org.Hs.eg.db)
#library(org.Ss.eg.db)
library(ggplot2)

inputdata <- read.table("Myocte_Findallmarkers.xls",header=T,sep="\t")
inputdata_filter <- inputdata[(inputdata$avg_log2FC > 1 & inputdata$p_val_adj < 0.05),]
head(inputdata_filter)


group <- data.frame(gene=inputdata$gene,group=inputdata$Broad)
Gene_ID <- bitr(geneID = inputdata$gene,fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")
anay_data <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')
data_GO <- compareCluster(ENTREZID~group, data=anay_data, fun="enrichGO", OrgDb="org.Hs.eg.db", ont = "BP",
                          pAdjustMethod = "BH",pvalueCutoff = 0.05,qvalueCutoff = 0.05)
write.table(data_GO,file="GOenrich.txt",sep = "\t")


pdf("dotplot.pdf",width = 8,height = 11)
dotplot(data_GO_sim, showCategory=5,font.size = 5)
dev.off()
