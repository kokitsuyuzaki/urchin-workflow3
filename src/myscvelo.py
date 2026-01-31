# -*- coding: utf-8 -*-
import sys
import scvelo as scv
import scanpy
import pandas as pd

args = sys.argv
mode = args[1]
infile1 = args[2]
infile2 = args[3]
outfile = args[4]

# Loading
adata1 = scanpy.read(infile1)
adata2 = scanpy.read_h5ad(infile2)

# Assign PCA, UMAP, and cluster ID calculated by Seurat
adata1.obs['clusters'] = pd.Categorical(adata2.obs['seurat_clusters'])
adata1.obsm['X_pca'] = adata2.obsm['X_pca']
adata1.obsm['X_umap'] = adata2.obsm['X_umap']
# 既存の座標・クラスタ移植
if "X_pca" in adata2.obsm:  
    adata1.obsm["X_pca"] = adata2.obsm["X_pca"].copy()
if "X_umap" in adata2.obsm: 
    adata1.obsm["X_umap"] = adata2.obsm["X_umap"].copy()

# クラスタを転送
if "seurat_clusters" in adata2.obs:
    adata1.obs["clusters"] = pd.Categorical(adata2.obs["seurat_clusters"])

# celltypeとcelltype_colorsを安全に転送
if "celltype" in adata2.obs:
    # 文字列に変換し、NAを処理
    celltype_values = adata2.obs["celltype"].values
    # numpy配列の場合の処理
    if hasattr(celltype_values, 'astype'):
        celltype_str = celltype_values.astype(str)
    else:
        celltype_str = [str(x) for x in celltype_values]
    adata1.obs["celltype"] = celltype_str

if "celltype_colors" in adata2.obs:
    # 同様に色情報も文字列化
    colors_values = adata2.obs["celltype_colors"].values
    if hasattr(colors_values, 'astype'):
        colors_str = colors_values.astype(str)
    else:
        colors_str = [str(x) for x in colors_values]
    adata1.obs["celltype_colors"] = colors_str

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
