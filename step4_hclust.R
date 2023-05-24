
library(Seurat)
library(tidyverse)


mg <- readRDS('merge_PigCellAtlas.RDS')
#*** select hvg 3000
hvgs_mat <- mg@assays$RNA@data[VariableFeatures(mg),]
write.table(row.names(hvgs_mat),file="hvgs_mat.txt",row.names=F,col.names=F,quote=F)
#*** CPM value
mg_counts <- mg@assays$RNA@counts
seuratCPM <- RelativeCounts(mg_counts, scale.factor = 1e6, verbose = TRUE)
seuratCPM_hvgs <- seuratCPM[row.names(hvgs_mat),]
dim(seuratCPM_hvgs)
log2_seuratCPM_hvgs <- log2(seuratCPM_hvgs+1)
dim(log2_seuratCPM_hvgs)
log2_seuratCPM_hvgs_mat <- as.matrix(log2_seuratCPM_hvgs)
log2_seuratCPM_hvgs_t <- t(log2_seuratCPM_hvgs_mat)
dim(log2_seuratCPM_hvgs_t)
log2_seuratCPM_hvgs_t_cluster <- aggregate(log2_seuratCPM_hvgs_t, by=list(mg@meta.data[["celltype"]]), FUN="mean") 
log2_seuratCPM_hvgs_t_cluster_t <- t(log2_seuratCPM_hvgs_t_cluster)
write.table(log2_seuratCPM_hvgs_t_cluster_t,file="log2_seuratCPM_hvgs3000.txt",sep="\t",quote=F,col.names = F)

#*** hclust
input<-read.table("log2_seuratCPM_hvgs3000.txt",sep="\t",header = T,row.names = 1)
input_t <- as.matrix(t(input))
hclust_data <- hclust(dist(input_t),method="ward.D2")
pdf("log2_seuratCPM_hvgs3000_ward.D2.pdf",width = 20,height = 10)
plot(hclust_data,hang=-1)
dev.off()

