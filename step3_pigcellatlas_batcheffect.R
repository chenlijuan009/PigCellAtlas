
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





