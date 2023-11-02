
library(Seurat)
library(ggplot2)
scRNA <- readRDS("Enterocyte.RDS")

ClusterColor = c("#e7a66b","#72f8d1","#771415","#c48c8b","#416f78","#ee87cd","#a0b8cb","#55a15c","#8247f5")
scRNA$celltype <- factor(scRNA$celltype,
                        levels = c("ISC_1","ISC_2","Enterocyte progenitor","Immat enterocyte_1",
                                   "Immat enterocyte_2", "Immat enterocyte_3",
                                   "Mat enterocyte","BEST4+ enterocyte","Enteroendocrine cell")) 

markergenes <- c("OLFM4","LGR5","MUC13","SI","FUT8","APOB","RAB3C","CHGA","LAMA3","SLC16A1","BEST4","STXBP5L")
markergenes <- unique(markergenes)
Idents(scRNA) <- scRNA@meta.data$celltype
pdf("Enetrocyte_markergene.pdf",width = 4,height=5)
VlnPlot(scRNA, features = markergenes,stack=TRUE,pt.size=0.5,fill.by='ident',flip=TRUE,cols=ClusterColor)+ 
       theme(axis.text.x=element_text(angle = 90), axis.ticks.x = element_line(),legend.position = 'none',
             axis.title.y = element_text(hjust = 0.5, angle = 0))+
       labs(x=NULL,y=NULL)
dev.off()


