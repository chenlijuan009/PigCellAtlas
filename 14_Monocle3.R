rm(list=ls())
library(monocle3)
library(dplyr)
library(ggplot2)
library(Seurat)
library(SeuratWrappers)

Enterocyte.data <-readRDS("Enterocyte.RDS")
Idents(Enterocyte.data)<-Enterocyte.data$celltype
#-----------------------------------------------------------------------------------------------------------------
#***Part1 mocnocle3 workflow
cds <- as.cell_data_set(Enterocyte.data)
cds
colData(cds) # to get cell metadata
fData(cds) # to gene metatdata
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name <- rownames(fData(cds))
counts(cds)
#-----------------------------------------------------------------------------------------------------------------
#***Part 2 pre-processing
cds <- preprocess_cds(cds,num_dim = 100)
plot_pc_variance_explained(cds)
cds <- reduce_dimension(cds)
plot_cells(cds)
cds <- cluster_cells(cds)
plot_cells(cds,color_cells_by = "partition")
list_cluster <- Enterocyte.data@active.ident
cds@clusters$UMAP$clusters <- list_cluster
cds@int_colData@listData$reducedDims$UMAP <- Enterocyte.data@reductions$umap@cell.embeddings

cluster.before.trajectory<-plot_cells(cds,
                                      color_cells_by = 'cluster',
                                      label_groups_by_cluster = FALSE,
                                      group_label_size = 3)+theme(legend.position = "right")
  
cluster_col=c("#e7a66b","#c48c8b","#416f78","#ee87cd","#55a15c","#a0b8cb","#8247f5")
cds$celltype <-factor(cds$celltype,
                      levels = c("ISC_1","Immat enterocyte_1","Immat enterocyte_2",
                                 "Immat enterocyte_3","BEST4+ enterocyte",
                                 "Mat enterocyte","Enteroendocrine cell")) 

cluster.names<-plot_cells(cds,color_cells_by = "celltype",label_groups_by_cluster = FALSE,group_label_size = 5)+
               scale_color_manual(values=cluster_col)+
               theme(legend.position = "none")
#-----------------------------------------------------------------------------------------------------------------
#***Part3 learn trajectory graph
library(tidydr)
cds <- learn_graph(cds,use_partition = FALSE)
p1  <- plot_cells(cds,color_cells_by='celltype',label_groups_by_cluster=FALSE,
                  label_branch_points = FALSE,label_roots = FALSE,label_leaves = FALSE,group_label_size = 5)+
       scale_color_manual(values=cluster_col)+
       theme_dr(xlength = 0.1,ylength = 0.1,arrow = grid::arrow(length = unit(0.1, "inches"),type = "closed"))+ 
       theme(aspect.ratio = 1,panel.grid = element_blank())+
       NoLegend()+
       theme(panel.background=element_blank(),panel.grid = element_blank(),
             plot.title = element_text(hjust=0.5,size=12),legend.text = element_text(size=12))
ggsave("plot_cells_cluster.pdf",p1,width=5,height=5)

#-----------------------------------------------------------------------------------------------------------------
#***Part4 order cells
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = colnames(cds[,clusters(cds)=="ISC_1"]))
p4  <- plot_cells(cds,color_cells_by = "pseudotime",
                  label_cell_groups = FALSE,
                  label_groups_by_cluster=FALSE,
                  label_leaves=FALSE,
                  label_branch_points=FALSE,
                  label_roots = FALSE)+ 
                 theme_dr(xlength = 0.1,ylength = 0.1,arrow = grid::arrow(length = unit(0.1, "inches"),type = "closed"))+ 
                 theme(aspect.ratio = 1,panel.grid = element_blank())+
                 theme(panel.background=element_blank(),panel.grid = element_blank(),
                       plot.title = element_text(hjust=0.5,size=12),legend.text = element_text(size=12))
ggsave("plot_cells_pseudotime.pdf",p4,width=5,height=5)     

#-----------------------------------------------------------------------------------------------------------------
#***Part5 cell ordered by monocel3 pesuodotime
pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
head(data.pseudo)
write.table(data.pseudo,file="data.pseudo.txt",sep="\t",quote=F)

p5<- ggplot(data.pseudo,aes(monocle3_pseudotime,reorder(celltype,monocle3_pseudotime,median),fill=celltype))+
     theme_classic()+
     geom_boxplot()+
     NoLegend()+
     scale_fill_manual(values=cluster_col)+
     theme(axis.text = element_text(colour = "black",size=12))+
     labs(x="",y="")
ggsave("plot_cells_pseudotime_box_median.pdf",p5,width=6,height=5)  

saveRDS(cds,file="monocle3_result.RDS")



