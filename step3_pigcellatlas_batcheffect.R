
library(Seurat)
library(reticulate) 
library(tidyverse)

use_miniconda("/.local/share/r-miniconda/bin/python")
scanorama <- import("scanorama")

Adipose <- readRDS("Adipose_decont_filter_singlet.so.rds")
Cerebellum <- readRDS("Cerebellum_decont_filter_singlet.so.rds")
Cerebrum <- readRDS("Cerebrum_decont_filter_singlet.so.rds")
Colon <- readRDS("Colon_decont_filter_singlet.so.rds")
Duodenum <- readRDS("Duodenum_decont_filter_singlet.so.rds")
Heart <- readRDS("Heart_decont_filter_singlet.so.rds")
Hypothalamus <- readRDS("Hypothalamus_decont_filter_singlet.so.rds")
Ileum <- readRDS("Ileum_decont_filter_singlet.so.rds")
Jejunum <- readRDS("Jejunum_decont_filter_singlet.so.rds")
Kidney <- readRDS("Kidney_decont_filter_singlet.so.rds")
Liver <- readRDS("Liver_decont_filter_singlet.so.rds")
Lymph <- readRDS("Lymph_decont_filter_singlet.so.rds")
Muscle <- readRDS("Muscle_decont_filter_singlet.so.rds")
Ovary <- readRDS("Ovary_decont_filter_singlet.so.rds")
Pancreas <- readRDS("Pancreas_decont_filter_singlet.so.rds")
Pituitary <- readRDS("Pituitary_decont_filter_singlet.so.rds")
Spleen <- readRDS("Spleen_decont_filter_singlet.so.rds")
Testis <- readRDS("Testis_decont_filter_singlet.so.rds")
Uterus <- readRDS("Uterus_decont_filter_singlet.so.rds")

mg<-merge(Adipose,y=c(Cerebellum,Cerebrum,Colon,Duodenum,Heart,Hypothalamus,Ileum,Jejunum,
                      Kidney,Liver,Lymph,Muscle,Ovary,Pancreas,Pituitary,Spleen,Testis,Uterus),
          add.cell.ids=c("Adipose","Cerebellum","Cerebrum","Colon","Duodenum","Heart","Hypothalamus",
                         "Ileum","Jejunum","Kidney","Liver","Lymph","Muscle","Ovary","Pancreas",
                         "Pituitary","Spleen","Testis","Uterus")
          )

sample <- data.frame(sample=c(rep("Adipose",ncol(Adipose)),rep("Cerebellum",ncol(Cerebellum)),
                              rep("Cerebrum",ncol(Cerebrum)),rep("Colon",ncol(Colon)),
                              rep("Duodenum",ncol(Duodenum)),rep("Heart",ncol(Heart)),
                              rep("Hypothalamus",ncol(Hypothalamus)),rep("Ileum",ncol(Ileum)),
                              rep("Jejunum",ncol(Jejunum)),rep("Kidney",ncol(Kidney)),
                              rep("Liver",ncol(Liver)),rep("Lymph",ncol(Lymph)),
                              rep("Muscle",ncol(Muscle)),rep("Ovary",ncol(Ovary)),
                              rep("Pancreas",ncol(Pancreas)),rep("Pituitary",ncol(Pituitary)),
                              rep("Spleen",ncol(Spleen)),rep("Testis",ncol(Testis)),
                              rep("Uterus",ncol(Uterus))
                              ))
table(sample)
rownames(sample) <- row.names(mg@meta.data)
mg<- AddMetaData(mg, sample)
str(mg)

mg[["percent.mt"]] <- PercentageFeatureSet(mg, pattern = "^MT-") 
mg$log10GenesPerUMI <- log10(mg$nFeature_RNA)/log10(mg$nCount_RNA)

#*** filter 
minGene=200
maxGene=5000
pctMT=5 #Remove cells with more than 5% mitochondrial reads.
minCount=500
maxCount=15000
max_log10GenesPerUMI=0.8
mg <- subset(mg, subset = nFeature_RNA > minGene & 
                          nFeature_RNA < maxGene &
                          log10GenesPerUMI > max_log10GenesPerUMI &
                          nCount_RNA> minCount & 
                          nCount_RNA<maxCount & 
                          percent.mt < pctMT
             )
mg <- NormalizeData(mg) %>% 
      FindVariableFeatures(nfeatures = 3000) %>%
      ScaleData() %>%
      RunPCA(npcs = 50, verbose=F) 
#*** Select PCA
pct <- mg[['pca']]@stdev / sum( mg[['pca']]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs 
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
# Elbow plot to visualize 
png("select_pca_value_elbow.png", width=9, height=6, res=300, units='in')
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
      geom_text() + 
      geom_vline(xintercept = 90, color = "grey") +
      geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
      theme_bw()
dev.off()
#*** Scanorama 
ifnb.list <- SplitObject(mg, split.by = "sample")
assaylist <- list()
genelist <- list()
for(i in 1:length(ifnb.list))
                            {
                             assaylist[[i]] <- t(as.matrix(GetAssayData(ifnb.list[[i]], "data")))
                             genelist[[i]] <- rownames(ifnb.list[[i]])
                            }
lapply(assaylist, dim)
mat_integrated <- scanorama$integrate(assaylist, genelist)
integrated.corrected.data <- scanorama$correct(assaylist, genelist,return_dimred=TRUE, return_dense=TRUE)
dr_scanorama <- do.call(rbind, integrated.corrected.data[[1]]) 
rownames(dr_scanorama) <- do.call(c, lapply(assaylist, rownames))
colnames(dr_scanorama) <- paste0("PC_", 1:100)
stdevs <- apply(dr_scanorama, MARGIN = 2, FUN = sd)
print(stdevs)
dr_scanorama <- dr_scanorama[colnames(mg),]
mg[['scanorama']] <- CreateDimReducObject(dr_scanorama, key="SCANORAMA_")
saveRDS(mg, "seurat_scanorma.rds")

mg <- FindNeighbors(mg, reduction='scanorama', dims = 1:pcs)
res.num <- seq(0.1,1.5,by=0.2)
plot_list_umap <- list()
plot_sample_umap <- list()
#*** iteration  UMAP
pdf("umapplot.pdf",width = 25,height = 13)
for (i in res.num)
                 {
                  mg <- FindClusters(mg, resolution = i)
                  mg <- RunUMAP(object=mg, reduction='scanorama',dims = 1:pcs) 
                  plot_list_umap[[i]] <- DimPlot(mg, reduction = "umap",label=T,pt.size = .2)
                  plot_sample_umap[[i]] <- DimPlot(mg, reduction = "umap",group.by = "sample", pt.size = .2)
                  print(plot_list_umap[[i]]+plot_sample_umap[[i]])
                 }
dev.off()
#*** iteration  tSNE
plot_list_tsne <-list()
plot_sample_tsne <-list()
pdf("tsneplot.pdf",width = 25,height = 13)
for (i in res.num)
                  {
                   mg<- FindClusters(mg, resolution = i)
                   mg <- RunTSNE(mg, reduction='scanorama', dims = 1:pcs) 
                   plot_list_tsne[[i]] <- DimPlot(mg, reduction = "tsne", label = TRUE, pt.size = .1)
                   plot_sample_tsne[[i]] <- DimPlot(mg, reduction = "tsne", group.by = "sample", pt.size = .1)
                   print(plot_list_tsne[[i]]+plot_sample_tsne[[i]])
                  }
dev.off()
saveRDS(mg, 'mg_scanorama_umap_tsne.rds')
#*** Find markers
Markers <- FindAllMarkers(mg, min.pct = 0.25, logfc.threshold = 0.25)
write.table(Markers,file="findallmarkers.xls",sep="\t")
saveRDS(Markers , 'findallmarkers_table.RDS')

options(repr.plot.height = 4, repr.plot.width = 6)
p3<-DimPlot(mg, reduction = "umap", group.by = "sample", pt.size = .01, split.by = 'sample')
ggsave("umap_split_sample.pdf", plot = p3, width = 30, height = 13,units = "cm") 
p4<-DimPlot(mg, reduction = "umap", group.by = "sample", pt.size = .01)
p5<-DimPlot(mg, reduction = "umap", label = TRUE, pt.size = .01)
ggsave("umap_sample_plot.pdf", plot = p4+p5, width = 30, height = 15,units = "cm") 
p6<-DimPlot(mg, reduction = "tsne", group.by = "sample", pt.size = .01)
p7<-DimPlot(mg, reduction = "tsne", label = TRUE, pt.size = .01)
ggsave("tsne_sample_plot.pdf", plot = p6+p7, width = 30, height = 15,units = "cm") 
p8<-DimPlot(mg, reduction = "tsne", group.by = "sample", pt.size = .1, split.by = 'sample')
ggsave("tsne_split_sample_plot.pdf", plot = p8, width = 30, height = 13,units = "cm") 



