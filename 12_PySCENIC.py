

libary(Seurat)
#--------------------------------------------------------------------------------------------------------------------
#*** step1 count matrix data
mg<-readRDS("/Users/lijuanchen/R/pig180D/Myocyte_Findallmaker.RDS")
write.csv(t(as.matrix(mg@assays$RNA@counts)),file = "Myocyte_RNAcounts.csv")
#--------------------------------------------------------------------------------------------------------------------
#***step 2 counts matrix 2 loom file
import numpy as np
import scanpy as sc
import loompy as lp
import pandas as pd

import os, sys
os.getcwd()
os.listdir(os.getcwd()) 
myocyte=sc.read_csv("Myocyte_RNAcounts.csv");
myocyte.shape

row_attrs = {"Gene": np.array(myocyte.var_names),}
col_attrs = {
             "CellID": np.array(myocyte.obs_names),
             "nGene": np.array(np.sum(myocyte.X.transpose()>0, axis=0)).flatten(), 
             "nUMI": np.array(np.sum(myocyte.X.transpose(),axis=0)).flatten(),
             }
lp.create("Myocyte_RNAcounts.loom",myocyte.X.transpose(),row_attrs, col_attrs)
#--------------------------------------------------------------------------------------------------------------------
#***step3 download TF and cisTarget data
#transcription factors
wget https://raw.githubusercontent.com/aertslab/pySCENIC/master/resources/hs_hgnc_tfs.txt
# motif to TF annotation database
wget https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
# genome ranking database
wget https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/ gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather

#--------------------------------------------------------------------------------------------------------------------
#***step4 run pyscenic

pyscenic grn Myocyte_RNAcounts.loom \
/public/agis/yiguoqiang_group/chenlijuan/database/pyscencic_DB/pySCENIC-0.11.2/resources/hs_hgnc_tfs.txt \
-o adj.tsv \
--method grnboost2 \
--num_workers 10

pyscenic ctx adj.tsv \
    /public/agis/yiguoqiang_group/chenlijuan/database/pyscencic_DB/cistarget/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
    --annotations_fname /public/agis/yiguoqiang_group/chenlijuan/database/pyscencic_DB/motif2TF/motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
    --expression_mtx_fname Myocyte_RNAcounts.loom \
    --mode "dask_multiprocessing" \
    --output reg.csv \
    --mask_dropouts \
    --num_workers 10


pyscenic aucell \
    Myocyte_RNAcounts.loom \
    reg.csv \
    --output pyscenic_output.loom \
    --num_workers 10


