#!/usr/local/bin/Rscript Rscript
library(Seurat)
library(cowplot)
library(harmony)
library(dplyr)
library(tidyverse)

com_mg_Human <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Human_PBMC_Lung.RDS")
com_mg_Monkey <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Monkey_PBMC_Lung.RDS")
com_mg_Mouse <- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Mouse_PBMC_Lung.RDS")
com_mg_PigOursudy <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_PigOurstudy.RDS")
com_mg_PigWang <-readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_Wangstudy.RDS")

com_mg_Human_Kidney <- subset(com_mg_Human,subset=tissue %in% c("Adult_Kidney2","Adult_Kidney3","Adult_Kidney4_1","Adult_Kidney4_2"))
com_mg_Monkey_Kidney <-subset(com_mg_Monkey,subset=tissue %in% c("Kidney"))
com_mg_Mouse_Kidney <-subset(com_mg_Mouse,subset=tissue %in% c("Kidney1","Kidney2"))
com_mg_PigOursudy_Kidney <-subset(com_mg_PigOursudy,subset=tissue %in% c("Kidney"))
com_mg_PigWang_Kidney <-subset(com_mg_PigWang,subset=tissue %in% c("Kidney1","Kidney2","Kidney3","Kidney4"))


com_cross_species_Kidney <- merge(com_mg_Human_Kidney,y=c(com_mg_Monkey_Kidney,com_mg_Mouse_Kidney,com_mg_PigOursudy_Kidney,com_mg_PigWang_Kidney))

#-----------------------------------------------------------------------------------------------------------------------------------------------
#Run Harmony

com_cross_species<- readRDS("/public/agis/yiguoqiang_group/chenlijuan/anlysis_project/SC_tussie/cross-species/com_mg_PBMC_Lung/mg_Cross_species.RDS")
table(com_cross_species$Species_Tissue)
com_cross_species_Kidney <- subset(com_cross_species,subset=Species_Tissue %in% c("Human_Adult_Kidney2","Human_Adult_Kidney3",
                                                                                  "Human_Adult_Kidney4_1","Human_Adult_Kidney4_2",
                                                                                  "Monkey_Kidney","Mouse_Kidney1","Mouse_Kidney2","PigOurstudy_Kidney",
                                                                                  "PigWangstudy_Kidney1","PigWangstudy_Kidney2",
                                                                                  "PigWangstudy_Kidney3","PigWangstudy_Kidney4"))
table(com_cross_species_Kidney$Species_Tissue)
saveRDS(com_cross_species_Kidney,file="com_cross_species_Kidney.RDS")


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

png("scRNA_harmony_10X_select_pca_value_elbow.png", width=9, height=6, res=300, units='in')
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + geom_text() + 
  geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +theme_bw()
dev.off()

scRNA_harmony <- scRNA_harmony %>%  RunHarmony("Species_Tissue", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(scRNA_harmony, 'harmony')
scRNA_harmony <- scRNA_harmony %>% RunTSNE(reduction = "harmony", dims = 1:pcs) %>% 
                 RunUMAP(reduction = "harmony", dims = 1:pcs) %>% 
                 FindNeighbors(reduction = "harmony", dims = 1:pcs) %>% 
                 FindClusters(resolution = c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4)) 

saveRDS(scRNA_harmony,"com_mg_cross_species_kidney_harmony.RDS")

