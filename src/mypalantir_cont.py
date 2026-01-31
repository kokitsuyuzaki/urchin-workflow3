# -*- coding: utf-8 -*-
import sys, os, pickle
import numpy as np
import pandas as pd
import scanpy as sc
import palantir

# args
infile  = sys.argv[1]
outfile = sys.argv[2]  # CSV想定

# ---- Load
adata = sc.read_h5ad(infile)

# ---- Preprocess（Palantir用に最小限）
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']].copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=50)

# ---- Diffusion maps → multiscale
dm_res = palantir.utils.run_diffusion_maps(
    pd.DataFrame(adata.obsm['X_pca'], index=adata.obs_names),
    n_components=30
)
ms_res = palantir.utils.determine_multiscale_space(dm_res)

# ---- 24h各細胞型から代表セルを選ぶ（PCA重心に最も近い細胞）
early_cells = []
if {"sample", "celltype"}.issubset(adata.obs.columns):
    m24 = adata.obs["sample"].astype(str) == "cont-24h"
    if m24.any():
        sub = adata.obs.loc[m24].copy()
        for ct, idx_cells in sub.groupby(sub["celltype"].astype(str)).groups.items():
            pos = adata.obs_names.get_indexer(pd.Index(idx_cells))
            X = adata.obsm["X_pca"][pos]
            ctr = X.mean(axis=0)
            i0 = np.argmin(((X - ctr) ** 2).sum(axis=1))
            early_cells.append(adata.obs_names[pos[i0]])

# フォールバック（条件に合うものが無い場合は先頭）
if not early_cells:
    early_cells = [adata.obs_names[0]]

# ---- Palantirを代表セルごとに実行 → 擬似時間を0-1正規化してマージ（中央値）
pts = []
for ec in early_cells:
    pr = palantir.core.run_palantir(ms_res, early_cell=ec)
    pt = pr.pseudotime.reindex(adata.obs_names)
    pt = (pt - pt.min()) / (pt.max() - pt.min())
    pts.append(pt.rename(ec))

pt_df = pd.concat(pts, axis=1)
pt_consensus = pt_df.median(axis=1).rename("palantir_pseudotime")

# ---- Save
pt_consensus.to_frame().to_csv(outfile, index=True)
