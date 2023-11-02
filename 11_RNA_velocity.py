#-----------------------------------------------------------------------------------------------------------
#***Step1 Export seurat data
library(Seurat)
library(Matrix)
output='/Users/lijuanchen/R/pig180D/Muscle/'
seurat_obj <-readRDS("/Users/lijuanchen/R/pig180D/Muscle/Muscle_abbr.RDS")
Idents(seurat_obj) <- seurat_obj@meta.data$celltype
seurat_obj$barcode <- colnames(seurat_obj)
seurat_obj$UMAP_1 <- seurat_obj@reductions$umap@cell.embeddings[,1]
seurat_obj$UMAP_2 <- seurat_obj@reductions$umap@cell.embeddings[,2]
seurat_obj$TSNE_1 <- seurat_obj@reductions$tsne@cell.embeddings[,1]
seurat_obj$TSNE_2 <- seurat_obj@reductions$tsne@cell.embeddings[,2]
write.csv(seurat_obj@meta.data, file=paste0(output, 'metadata.csv'), quote=F, row.names=F)
counts_matrix <- GetAssayData(seurat_obj, assay='RNA', slot='counts')
writeMM(counts_matrix, file=paste0(output, 'counts.mtx'))
write.csv(seurat_obj@reductions$pca@cell.embeddings,file=paste0(output, 'pca.csv'), quote=F, row.names=F) 
write.table(data.frame('gene'=rownames(counts_matrix)),file=paste0(output, 'gene_names.csv'),quote=F,row.names=F,col.names=F)            
#-----------------------------------------------------------------------------------------------------------
#***Step3 convert h5ad data
import scanpy as sc
import anndata
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix
import numpy as np
import os
import pandas as pd
os.chdir("/Users/lijuanchen/R/pig180D/Muscle/RNAveocity/")
os.getcwd()
# load sparse matrix:
X = io.mmread("counts.mtx")
# create anndata object
adata = anndata.AnnData(
    X=X.transpose().tocsr()
)
# load cell metadata:
cell_meta = pd.read_csv("metadata.csv")
# load gene names:
with open("gene_names.csv", 'r') as f:
    gene_names = f.read().splitlines()
# set anndata observations and index obs by barcodes, var by gene names
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']
adata.var.index = gene_names
# load dimensional reduction:
pca = pd.read_csv("pca.csv")
pca.index = adata.obs.index

# set pca and umap
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'].to_numpy(), adata.obs['UMAP_2'].to_numpy())).T



#-----------------------------------------------------------------------------------------------------------
#***Step4 Merge velocity data

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=300, frameon=False)
cr.settings.verbosity = 2
ldata1 = scv.read('/Users/lijuanchen/R/pig180D/Muscle/RNAveocity/possorted_genome_bam_TR7I8.loom', cache=True)
barcodes = [bc.split(':')[1] for bc in ldata1.obs.index.tolist()]
barcodes = [bc[0:len(bc)-1] + '_Muscle' for bc in barcodes]
ldata1.obs.index = barcodes
ldata1.var_names_make_unique()
adata = scv.utils.merge(adata, ldata1)
adata.write('muscle_data.h5ad')


#-----------------------------------------------------------------------------------------------------------
#***Step5 run RNA_velocity
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import scipy as sci
import scipy.sparse
import seaborn as sns
import pandas as pd
import numpy as np
import itertools
from collections import Counter
import re
import glob
import os
import scvelo as scv

scv.settings.verbosity = 3  # show errors(0), warnings(1), info(2), hints(3)
scv.settings.presenter_view = True  # set max width size for presenter view
scv.set_figure_params('scvelo')  # for beautified visualization
import scanpy as sc
import anndata
from anndata import read_h5ad
from anndata import AnnData
from scipy import io
from scipy.sparse import coo_matrix, csr_matrix

#os.chdir("/Users/lijuanchen/R/pig180D/resultfigure/RNAvelocity_Muscle")
os.getcwd()
adata= sc.read_h5ad('muscle_data.h5ad')

# Preprocess the data
scv.pl.proportions(adata)
scv.pp.filter_genes(adata, min_shared_counts=20)
scv.pp.normalize_per_cell(adata)
scv.pp.filter_genes_dispersion(adata, n_top_genes=2000)
scv.pp.log1p(adata)
scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata, svd_solver='arpack')
sc.pp.neighbors(adata,n_neighbors=30)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
# Estimate RNA velocity

scv.tl.recover_dynamics(adata)
scv.tl.velocity(adata, mode='dynamical')
scv.tl.velocity_graph(adata)
scv.settings.set_figure_params('scvelo', dpi=300, dpi_save=300)
scv.pl.velocity_embedding_stream(adata, basis='umap',color ='celltype', legend_loc='right',save='velocity_embedding_stream.pdf')
scv.settings.set_figure_params('scvelo', dpi=300, dpi_save=300)
scv.pl.velocity_embedding_grid(adata, basis='umap', color ="celltype", legend_loc='right',save='velocity_embedding_grid.pdf')

scv.tl.velocity_confidence(adata)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata, c=keys, cmap='coolwarm', perc=[5, 95])
df = adata.obs.groupby('celltype')[keys].mean().T
df.style.background_gradient(cmap='coolwarm', axis=1)
scv.pl.scatter(adata, color='velocity_confidence', cmap='coolwarm', perc=[5, 95],save="velocity_confidence.pdf")

#Identify important genes
scv.tl.rank_velocity_genes(adata, groupby='celltype', min_corr=.3)

df = scv.DataFrame(adata.uns['rank_velocity_genes']['names'])
df.head()
important=df.iloc[0,:].tolist()
important = important + df.iloc[1,:].tolist() + df.iloc[2,:].tolist()
scv.settings.set_figure_params('scvelo', dpi=300, dpi_save=300)
scv.pl.velocity(adata, important,save="velocity_important.pdf")

# Latent time
scv.tl.latent_time(adata)
scv.pl.scatter(adata, color='latent_time', color_map='gnuplot', size=80,save="scatter_latenttime.pdf")


# heatmap 
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adata, var_names=top_genes, sortby='latent_time', col_color='celltype', n_convolve=100,save="latent_time_heatmap.pdf")

# top-likelihood genes
top_genes = adata.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adata, basis=top_genes[:15], ncols=5, frameon=False,save="top_genes_scatter.pdf")

# cluster-specific top-likelihood gens
scv.tl.rank_dynamical_genes(adata, groupby='celltype')
df = scv.get_df(adata, 'rank_dynamical_genes/names')
df.head(5)
adata.write("muscle_data_result.h5ad")


