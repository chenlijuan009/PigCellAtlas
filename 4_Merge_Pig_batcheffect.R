#!/usr/local/bin/Rscript Rscript
library(Seurat)
library(tidyverse)
#-----------------------------------------------------------------------------------------------------------------------------------------------
#*** step1 Merge pig Ourstudy data
setwd("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/")
listRDS <- dir(pattern = "*.rds") 
print(length(listRDS))

object_list <- lapply(X = listRDS, FUN = function(x) { x <- readRDS(x) })
mg <- merge(object_list[[1]], y = object_list[2:length(object_list)],
                        add.cell.ids=c("Adipose","Cerebellum","Cerebrum","Colon","Duodenum","Heart","Hypothalamus","Ileum","Jejunum",
                                       "Kidney","Liver","Lymph","Muscle","Ovary","Pancreas","Pituitary","Spleen","Testis","Uterus")
                            

species <- data.frame(species =c(rep('PigOurstudy',dim(mg@meta.data)[1])))
rownames(species) <- row.names(mg@meta.data)
mg<- AddMetaData(mg,species)

saveRDS(mg, 'mg_PigOurstudy.RDS')
#-----------------------------------------------------------------------------------------------------------------------------------------------
#*** step2 step1 Merge pig Wangstudy data

listRDS <- dir(pattern = "*.rds") 
print(length(listRDS))

object_list <- lapply(X = listRDS, FUN = function(x) { x <- readRDS(x) })


mgmgWang <- merge(object_list[[1]], y = object_list[2:length(object_list)],
                        add.cell.ids=c("AdiposeS","AdiposeV","AreaPostrema1","AreaPostrema2","AreaPostrema3","Brain","Cerebellum1","Cerebellum2","Cerebellum3",
                                       "Cerebellum4","Cerebellum5","Cerebellum6","Cerebellum7","Cerebellum8","Heart1","Heart2","Heart3","Heart4","Heart5","Intestine",
                                       "Kidney1","Kidney2","Kidney3","Kidney4","Liver1","Liver2","Liver3","Liver4"，"LiverW","Lung","OVoLT1","OVoLT2","OVoLT3","OVoLT4",
                                       "PBMC","Retina","Retina1","Retina2","Spleen1","Spleen2","Spleen3","SpleenW","SubfomicalOrgan1","SubfomicalOrgan2","SubfomicalOrgan3"))


species <- data.frame(species =c(rep('PigWangstudy',dim(mg@meta.data)[1])))
rownames(sample) <- row.names(mgWang@meta.data)
mgWang<- AddMetaData(mgWang, sample)
str(mgWang)
rownames(species) <- row.names(mgWang@meta.data)
mgWang<- AddMetaData(mgWang, species)

saveRDS(mgWang, 'mg_Wangstudy.RDS')

mgourWang <- merge(mg,mgWang,add.cell.ids=c("Pig_Ourstudy","Pig_Wangstudy"))
saveRDS(mgourWang,file="mg_ourwangstudy.RDS")
#-------------------------------------------------------------------------------------------------------------
#***Step3 run harmony
scRNA_harmony <- readRDS("mg_OurWangstudy.RDS")
minGene=200
maxGene=5000
pctMT=30 #Remove cells with more than 5% mitochondrial reads.
# pctMT=5
minCount=500
maxCount=15000
max_log10GenesPerUMI=0.8
scRNA_harmony <- subset(scRNA_harmony, subset = nFeature_RNA > minGene & nFeature_RNA < maxGene & log10GenesPerUMI > max_log10GenesPerUMI 
                        & nCount_RNA> minCount & nCount_RNA<maxCount & percent.mt < pctMT)

scRNA_harmony <- scRNA_harmony %>% NormalizeData(verbose = FALSE) %>%
                 FindVariableFeatures(selection.method = "vst", nfeatures = 3000) %>% 
                 ScaleData(verbose = FALSE) %>% 
                 RunPCA(pc.genes = scRNA_harmony@var.genes, npcs = 100, verbose = FALSE)

pct <- scRNA_harmony[['pca']]@stdev / sum( scRNA_harmony[['pca']]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs 
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
# Elbow plot to visualize 
png("scRNA_harmony_10X_select_pca_value_elbow.png", width=9, height=6, res=300, units='in')
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + geom_text() + 
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +theme_bw()
dev.off()

scRNA_harmony <- scRNA_harmony %>%  RunHarmony("sample", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(scRNA_harmony, 'harmony')
scRNA_harmony <- scRNA_harmony %>%
                 RunTSNE(reduction = "harmony", dims = 1:pcs) %>% 
                 RunUMAP(reduction = "harmony", dims = 1:pcs) %>% 
                 FindNeighbors(reduction = "harmony", dims = 1:pcs) %>% 
                 FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) 
saveRDS(scRNA_harmony,"mg_OurWangstudy_harmony.RDS")


