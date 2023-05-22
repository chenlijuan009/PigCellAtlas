library(magrittr)
library(Seurat)
library(celda)
library(scater)
library(DoubletFinder)
#***Ambient RNA removal
Adipose <- Read10X(data.dir='/Adipose/filtered_feature_bc_matrix')
Adipose.sce <- SingleCellExperiment(list(counts = Adipose))
Adipose.sce <- decontX(Adipose.sce)
Adipose_decont.sce <- Adipose.sce[,which(Adipose.sce$decontX_contamination < 0.8)]
Adipose_decont.so <- CreateSeuratObject(round(decontXcounts(Adipose_decont.sce)))
Adipose_obj = Adipose_decont.so
#*** decontXcounts
Adipose_obj[['percent.mt']] <- PercentageFeatureSet(object = Adipose_obj, pattern = '^MT-')
fivenum(Adipose_obj[['percent.mt']][,1])
Adipose_obj@meta.data$log10GenesPerUMI <- log10(Adipose_obj@meta.data$nFeature_RNA)/log10(Adipose_obj@meta.data$nCount_RNA)
saveRDS(Adipose_obj,'Adipose_decont_qc.so.rds')

#***filter
#*** minGene=200 & maxGene=5000 & pctMT=5 & minCount=500 & maxCount=15000
Adipose_obj <- subset(x = Adipose_obj, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & nCount_RNA > 500 & nCount_RNA < 15000 & percent.mt < 5 & log10GenesPerUMI > 0.8)
#*** scacle 
Adipose_obj <- NormalizeData(object = Adipose_obj, normalization.method = 'LogNormalize', scale.factor = 10000)
Adipose_obj <- FindVariableFeatures(object=Adipose_obj,selection.method = 'vst',nfeatures = 3000)
Adipose_obj_scal<- ScaleData(object =Adipose_obj,features = VariableFeatures(Adipose_obj))
Adipose_obj_scal<- RunPCA(object = Adipose_obj_scal, npcs = 50, pc.genes = VariableFeatures(object = Adipose_obj_scal))
#***Select number of pca
pct <- Adipose_obj_scal [['pca']]@stdev / sum(Adipose_obj_scal [['pca']]@stdev) * 100
cumu <- cumsum(pct)
co1 <- which(cumu > 90 & pct < 5)[1]
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
pcs <- min(co1, co2)
pcs
plot_df <- data.frame(pct = pct, cumu = cumu, rank = 1:length(pct))
png("select_pca_value_elbow.png", width=9, height=6, res=300, units='in')
ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
  geom_text() + geom_vline(xintercept = 90, color = "grey") +
  geom_hline(yintercept = min(pct[pct > 5]), color = "grey") +
  theme_bw()
dev.off()
#***UMAP
Adipose_obj_scal <- RunUMAP(Adipose_obj_scal, dims = 1:pcs)
#**TSNE
Adipose_obj_scal <- RunTSNE(Adipose_obj_scal, dims = 1:pcs)
saveRDS(Adipose_obj_scal,'Adipose_decont_filter_scal.so.rds')

#***detect doublet
sweep.res.list <- paramSweep_v3(Adipose_obj_scal, PCs = 1:pcs, sct = F)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pK_bcmvn <- bcmvn$pK[which.max(bcmvn$BCmetric)] %>% as.character() %>% as.numeric()
annotations <- Adipose_obj_scal@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)
DoubletRate = 0.075
nExp_poi <- round(DoubletRate*nrow(Adipose_obj_scal@meta.data))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
Adipose_doubletfinder_obj <- doubletFinder_v3(Adipose_doublet_obj, PCs = 1:pcs, pN = 0.25, pK = pK_bcmvn, nExp = nExp_poi.adj, reuse.pANN = paste0('pANN_0.25_',pK_bcmvn,'_',nExp_poi.adj), sct = F)
Adipose_decont_filter_singlet <- subset(Adipose_doubletfinder_obj, cells= rownames(Adipose_doubletfinder_obj@meta.data[Adipose_doubletfinder_obj@meta.data[[7]]=='Singlet',]))
saveRDS(Adipose_decont_filter_singlet, 'Adipose_decont_filter_singlet.so.rds')







































