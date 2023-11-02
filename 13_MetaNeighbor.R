library(MetaNeighbor)
library(scRNAseq)
library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(stringr)
#--------------------------------------------------------------------------------------------------------
#***step1   Seurat convert to SingleCellExpreiment

cal_pseudo_cell <- readRDS("cal_pseudo_cell.RDS")
cal_pseudo_cell_sce <- SingleCellExperiment(assays = list(counts=cal_pseudo_cell@assays$RNA@counts),
                                            colData=cal_pseudo_cell@meta.data)
colnames(colData(cal_pseudo_cell_sce))

study_id <- cal_pseudo_cell$species
celltype <- cal_pseudo_cell$celltype

var_genes = variableGenes(dat = cal_pseudo_cell_sce, exp_labels = study_id )

length(var_genes)
AUROC_scores = MetaNeighborUS(var_genes = var_genes,dat = cal_pseudo_cell_sce,
                              study_id = study_id ,cell_type = celltype)
head(AUROC_scores)  
write.table(AUROC_scores,file="AUROC_scores.txt",sep="\t",quote=F)
#fillter threshold
top_hits = topHits(cell_NV = AUROC_scores,dat = cal_pseudo_cell_sce,study_id = cal_pseudo_cell_sce$study_id,
                   cell_type = cal_pseudo_cell_sce$celltype,threshold = 0.9)
top_hits
write.table(top_hits,file="top_hits.txt",sep="\t",quote = F)

#--------------------------------------------------------------------------------------------------------
#***step2 draw heatmap
#***heatmap add labels

library(pheatmap)
AUCscore <- read.table("AUROC_scores.txt",sep = "\t",header=T,row.names = 1)
head(AUCscore)
celltype=read.table("celltype.txt",sep = "\t",header=T)
head(celltype)
Celltype <- celltype$Group.1

row_anno = data.frame(CellType= factor(Celltype),
                      species=factor(c(rep("Human",9),rep("Monkey",9),rep("Mouse",9),rep("Pig",9)))
                      )
col_anno = data.frame(CellType= factor(Celltype),
                      species=factor(c(rep("Human",9),rep("Monkey",9),rep("Mouse",9),rep("Pig",9)))
)
rownames(row_anno) <- rownames(AUCscore)                               
rownames(col_anno) <- colnames(AUCscore)

ann_colors = list(
                  species = c(Human = "#af33b2", Monkey = "#55a15c", Mouse = "#eae757",Pig="peachpuff"),
                  CellType= c(Endocrine="#af33b2",Endothelial="#55a15c",Epithelial="#eae757",Germline="#38659d",
                              Immune="#cf3b31",Islet="#e2c4ca",Muscle="#771415",Neural="#8247f5",Stromal="#9c7eaf")
)

pdf("pheatmap_heatmap.pdf",width =8,height = 7 )
pheatmap(AUCscore,scale = "none",
        clustering_method = "ward.D2", #euclidean ,ward.D 
        cluster_rows = TRUE,cluster_cols = TRUE,
        annotation_row = row_anno,
        annotation_col = col_anno,
        annotation_colors = ann_colors
        )
dev.off()



