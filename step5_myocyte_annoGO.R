library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Ss.eg.db)
library(ggplot2)
#*** each cell type
inputdata <- read.table("myocte_Findallmarkers_anno_broad_celltype_pos.xls",header=T,sep="\t")
inputdata_filter <- inputdata[(inputdata$avg_log2FC > 1 & inputdata$p_val_adj < 0.05),]
cell_cluster <- c('Muscle','SMC','CMC')
for( i in cell_cluster)
                       {
                         ls <- inputdata_filter %>% filter(cluster==i) 
                         write.table(ls ,paste('cell_cluster_',i,".txt",sep = ""),sep='\t',quote = F,row.names = F)
                         }
#*** Muscle
en_sym<-select(org.Hs.eg.db,keys = keys(org.Hs.eg.db),columns = c('ENSEMBL','SYMBOL','ENTREZID'))
Muscle <- read.table("cell_cluster_Muscle.txt",header=T,sep="\t")
com_Muscle <- en_sym[match(Muscle$gene,en_sym$SYMBOL),]

#**SMC 
SMC <- read.table("cell_cluster_SMC.txt",header=T,sep="\t")
com_SMC <- en_sym[match(SMC$gene,en_sym$SYMBOL),]

#**CMC
CMC <- read.table("cell_cluster_CMC.txt",header=T,sep="\t")
com_CMC <- en_sym[match(CMC$gene,en_sym$SYMBOL),]


com_Muscle_enterzID <- com_Muscle$ENTREZID[!is.na(com_Muscle$ENTREZID)]
com_SMC_enterzID <- com_SMC$ENTREZID[!is.na(com_SMC$ENTREZID)]
com_CMC_enterzID <- com_CMC$ENTREZID[!is.na(com_CMC$ENTREZID)]

#*** list cell type
gc<- list(Muscle = com_Muscle_enterzID, SMC = com_SMC_enterzID, CMC = com_CMC_enterzID)
gc
#*** enrichGO
celltype_enrichGO <- compareCluster(gc,fun="enrichGO",OrgDb="org.Hs.eg.db",ont= "BP",pvalueCutoff=1)
write.table(celltype_enrichGO,file="enrich_GO.txt",sep="\t",row.names = F)

enrich_dot <- dotplot(celltype_enrichGO, showCategory=5, includeAll=FALSE)
ggsave(plot = enrich_dot,filename = "enrich_dot.pdf",width = 8,height = 7)

#*** extract top 5 term each cluster
library(dplyr)
enrichGO <- read.table("enrich_GO.txt",sep="\t",header=T)
each_cluster <- enrichGO %>% group_by(Cluster) %>%  do(head(.,n=5))
write.table(each_cluster,file="enrich_each_cluster_GO.txt",sep="\t",row.names = F)











