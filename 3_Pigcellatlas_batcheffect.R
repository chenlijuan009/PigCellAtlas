
library(Seurat)
library(reticulate) 
library(tidyverse)

use_miniconda("/.local/share/r-miniconda/bin/python")
scanorama <- import("scanorama")
listRDS <- dir(pattern = "*.rds") 
print(length(listRDS))

object_list <- lapply(X = listRDS, FUN = function(x) { x <- readRDS(x) })


mg <- merge(object_list[[1]], y = object_list[2:length(object_list)],
          add.cell.ids=c("Adipose","Cerebellum","Cerebrum","Colon","Duodenum","Heart","Hypothalamus",
                         "Ileum","Jejunum","Kidney","Liver","Lymph","Muscle","Ovary","Pancreas",
                         "Pituitary","Spleen","Testis","Uterus")
          )

sample <- data.frame(sample =c(rep('Adipose',ncol(object_list[[1]])),rep('Cerebellum',ncol(object_list[[2]])),rep('Cerebrum',ncol(object_list[[3]])),
                               rep('Colon',ncol(object_list[[4]])),rep('Duodenum',ncol(object_list[[5]])),rep('Heart',ncol(object_list[[6]])),
                               rep('Hypothalamus',ncol(object_list[[7]])),rep('Ileum',ncol(object_list[[8]])),rep('Jejunum',ncol(object_list[[9]])),
                               rep('Kidney',ncol(object_list[[10]])),rep('Liver',ncol(object_list[[11]])),rep('Lymph',ncol(object_list[[12]])),
                               rep('Muscle',ncol(object_list[[13]])),rep('Ovary',ncol(object_list[[14]])),
                               rep('Pancreas',ncol(object_list[[15]])),rep('Pituitary',ncol(object_list[[16]])),
                               rep('Spleen',ncol(object_list[[17]])),rep('Testis',ncol(object_list[[18]])),rep('Uterus',ncol(object_list[[19]]))))

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
mg <- mg %>% NormalizeData() %>% 
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





