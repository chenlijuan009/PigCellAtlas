#!/usr/local/bin/Rscript Rscript
#-----------------------------------------------------------------------------
#***Part1 Merge HCL data
rm(list=ls())
library(Seurat)
library(tidyverse)

listRDS <- dir(pattern = "*.rds") 
print(length(listRDS))



mgHuman_object_list <- lapply(X = listRDS, FUN = function(x) { x <- readRDS(x) })

mgHuman <- merge(mgHuman_object_list[[1]], y = mgHuman_object_list[2:length(mgHuman_object_list)]))

species <- data.frame(species =c(rep('Human',dim(mg@meta.data)[1])))

rownames(species) <- row.names(mgHuman@meta.data)
mgHuman<- AddMetaData(mgHuman, species)


saveRDS(mgHuman, 'mg_Human.RDS')
                    

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#***Part2 Merge Monkey data
library(Seurat)
library(tidyverse)

listRDS <- dir(pattern = "*.rds") 
print(length(listRDS))

mgMonkey_object_list <- lapply(X = listRDS, FUN = function(x) { x <- readRDS(x) })

mgMonkey<-merge(mgHuman_object_list[[1]], y = mgHuman_object_list[2:length(mgHuman_object_list)]))

species <- data.frame(species =c(rep('Human',dim(mg@meta.data)[1])))

saveRDS(mgMonkey, 'mg_Monkey.RDS')



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#***Part3 Merge Mouse data
listRDS <- dir(pattern = "*.rds") 
print(length(listRDS))

mgMonkey_object_list <- lapply(X = listRDS, FUN = function(x) { x <- readRDS(x) })

mgMouse<-merge(mgHuman_object_list[[1]], y = mgHuman_object_list[2:length(mgHuman_object_list)]))

saveRDS(mgMouse, 'mg_Mous.RDS')


#----------------------------------------------------------------------------------------------------------------------------
#***Part3 Common_gene

GeneID <- read.table("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/gene_ID/final_comID.txt",header=T)
common_gene <- GeneID$Gene
mg_Human <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/mg_Human.RDS")
com_mg_Human <-mg_Human[common_gene,]
str(com_mg_Human)
saveRDS(com_mg_Human,file="com_mg_Human.RDS")

mg_Monkey <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/mg_Monkey.RDS")
com_mg_Monkey <- mg_Monkey[common_gene,]
str(com_mg_Monkey)
saveRDS(com_mg_Monkey,file="com_mg_Monkey.RDS")


mg_Mouse <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/mg_Mouse.RDS")
com_mg_Mouse <-mg_Mouse[common_gene,]
str(com_mg_Mouse)
saveRDS(com_mg_Mouse,file="com_mg_Mouse.RDS")

mg_PigOursudy <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/mg_PigOurstudy.RDS")
com_mg_PigOursudy <- mg_PigOursudy[common_gene,]
str(com_mg_PigOursudy)
saveRDS(com_mg_PigOursudy,file="com_mg_PigOurstudy.RDS")

mg_PigWang <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/mg_Wangstudy.RDS")
com_mg_PigWang <- mg_PigWang[common_gene,]
str(com_mg_PigWang)
saveRDS(com_mg_PigWang,file="com_mg_Wangstudy.RDS")

#----------------------------------------------------------------------------------------------------------------------------
#***Part4 Run harmony
com_mg_Human <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Human.RDS")
com_mg_Monkey <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Monkey.RDS")
com_mg_Mouse <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Mouse.RDS")
com_mg_PigOursudy <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_PigOursudy.RDS")
com_mg_PigWang <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_PigWang.RDS")

mg_Cross_species <- merge(com_mg_Human,y=c(com_mg_Monkey,com_mg_Mouse,com_mg_PigOursudy,com_mg_PigWang),add.cell.ids=c("Human","Monkey","Mouse","Pig","Pig"))

saveRDS(mg_Cross_species,file="mg_Cross_species.RDS")


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

scRNA_harmony <- scRNA_harmony %>%  RunHarmony("tissue", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(scRNA_harmony, 'harmony')
scRNA_harmony <- scRNA_harmony %>%
                 RunTSNE(reduction = "harmony", dims = 1:pcs) %>% 
                 RunUMAP(reduction = "harmony", dims = 1:pcs) %>% 
                 FindNeighbors(reduction = "harmony", dims = 1:pcs) %>% 
                 FindClusters(resolution = c(0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2, 2.2)) 

saveRDS(scRNA_harmony,"mg_Cross_species_harmony.RDS")

