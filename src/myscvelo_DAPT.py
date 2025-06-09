# -*- coding: utf-8 -*-
import sys
import scvelo as scv
import scanpy
import pandas as pd
import re

args = sys.argv

mode = args[1]
infile1 = args[2]
infile2 = args[3]
outfile = args[4]
# mode = "stochastic"
# infile1 = 'output/hpbase/aggr/velocyto/aggr.loom'
# infile2 = 'output/hpbase/integrated/seurat.h5ad'

# Loading
adata1 = scanpy.read(infile1)
adata2 = scanpy.read_h5ad(infile2)

# Filtering cells
cell_idx = adata2.obs.index[adata2.obs['sample'].str.contains(r"^DAPT-", regex=True)]
cell_pos = [adata2.obs_names.get_loc(idx) for idx in cell_idx]
adata1 = adata1[cell_pos].copy()
adata2 = adata2[cell_pos].copy()

# Assign PCA, UMAP, and cluster ID calculated by Seurat
adata1.obs['clusters'] = pd.Categorical(adata2.obs['seurat_clusters'])
adata1.obsm['X_pca'] = adata2.obsm['X_pca']
adata1.obsm['X_umap'] = adata2.obsm['X_umap']

# Preprocessing
scv.pp.filter_genes(adata1, min_shared_counts=20)
scv.pp.normalize_per_cell(adata1)
scv.pp.filter_genes_dispersion(adata1, n_top_genes=2000)
scv.pp.log1p(adata1)
scv.pp.filter_and_normalize(adata1, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata1, n_pcs=30, n_neighbors=1000)

# RNA Velocity
scv.tl.velocity(adata1, mode=mode)

# Velocity Graph
scv.tl.velocity_graph(adata1)

# Latent Time
if mode == "dynamical":
    scv.tl.recover_dynamics(adata1)
    scv.tl.latent_time(adata1)

# Save
adata1.write(outfile)
