# -*- coding: utf-8 -*-
# usage: python myscvelo_DAPT.py <stochastic|dynamical> <in_loom/h5ad> <ref_h5ad> <out_h5ad> [n_jobs]

import os, sys, re
import numpy as np, pandas as pd
import scanpy as sc, scvelo as scv
from scipy import sparse as sp

if len(sys.argv) < 5:
    raise SystemExit("usage: python myscvelo_DAPT.py <stochastic|dynamical> <in_loom/h5ad> <ref_h5ad> <out_h5ad> [n_jobs]")

mode, infile1, infile2, outfile = sys.argv[1:5]
n_jobs = int(sys.argv[5]) if len(sys.argv) > 5 else 8

n_pcs, n_neighbors, n_top_genes = 30, 30, 2000
sample_regex = r"^DAPT-"

os.environ["OMP_NUM_THREADS"]="1"
os.environ["OPENBLAS_NUM_THREADS"]="1"
os.environ["MKL_NUM_THREADS"]="1"
os.environ["NUMBA_NUM_THREADS"]=str(n_jobs)

adata1 = sc.read(infile1)
adata2 = sc.read_h5ad(infile2)

def base16(ix: pd.Index) -> pd.Index:
    s = pd.Index(ix.astype(str))
    b = s.str.extract(r'([ACGTN]{16})', expand=False)
    return pd.Index(b.fillna(""))

def loose_norm(ix: pd.Index) -> pd.Index:
    s = pd.Index(ix.astype(str))
    s = s.str.replace(r'^.*?:', '', regex=True)    # 接頭辞 "DAPT-24h:" を除去
    s = s.str.replace(r'[:].*$', '', regex=True)   # 残りの ":" 以降を除去
    s = s.str.replace(r'[_-]\d+$', '', regex=True) # "-1" や "_3" を除去
    s = s.str.replace(r'x$', '', regex=True)       # 末尾 'x' を除去（velocyto）
    return s

# Seurat側のサンプル抽出（無ければ全セル）
if "sample" in adata2.obs:
    keep = adata2.obs.index[adata2.obs["sample"].astype(str).str.contains(sample_regex, case=False, na=False)]
    if len(keep) == 0:
        keep = adata2.obs.index
else:
    keep = adata2.obs.index

# まず base16 で照合、ダメなら loose_norm
b1 = base16(adata1.obs_names)
b2 = base16(adata2.obs_names)
keep_b2 = base16(pd.Index(keep))
common = pd.Index(np.intersect1d(b1.values, keep_b2.values))

if len(common) == 0:
    l1 = loose_norm(adata1.obs_names)
    l2 = loose_norm(adata2.obs_names)
    keep_l2 = loose_norm(pd.Index(keep))
    common = pd.Index(np.intersect1d(l1.values, keep_l2.values))

if len(common) == 0:
    print("DEBUG adata1[:5]:", list(map(str, adata1.obs_names[:5])))
    print("DEBUG adata2[:5]:", list(map(str, adata2.obs_names[:5])))
    raise RuntimeError("共通セル0件。命名とフィルタを確認してください。")

# --- ここが修正点：重複キーを先勝ちで落としてから combine_first ---
m1_base = pd.Series(adata1.obs_names.values, index=b1)
m1_base = m1_base[~m1_base.index.duplicated(keep="first")]
m1_loose = pd.Series(adata1.obs_names.values, index=loose_norm(adata1.obs_names))
m1_loose = m1_loose[~m1_loose.index.duplicated(keep="first")]
map1 = m1_base.combine_first(m1_loose)

m2_base = pd.Series(adata2.obs_names.values, index=b2)
m2_base = m2_base[~m2_base.index.duplicated(keep="first")]
m2_loose = pd.Series(adata2.obs_names.values, index=loose_norm(adata2.obs_names))
m2_loose = m2_loose[~m2_loose.index.duplicated(keep="first")]
map2 = m2_base.combine_first(m2_loose)

# 順序（重複除去）
order = pd.Index(pd.unique(common))
adata1 = adata1[map1.reindex(order).values].copy()
adata2 = adata2[map2.reindex(order).values].copy()
adata1.obs_names_make_unique(); adata2.obs_names_make_unique()

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

# spliced/unspliced を疎に
for k in ("spliced","unspliced","ambiguous"):
    if k in adata1.layers and adata1.layers[k] is not None and not sp.issparse(adata1.layers[k]):
        adata1.layers[k] = sp.csr_matrix(adata1.layers[k])

# 前処理（重複適用しない）
scv.pp.filter_and_normalize(adata1, min_shared_counts=20, n_top_genes=n_top_genes, enforce=False, log=False)
scv.pp.moments(adata1, n_pcs=n_pcs, n_neighbors=n_neighbors)

# 速度
if mode == "dynamical":
    scv.tl.recover_dynamics(adata1, n_jobs=n_jobs, max_iter=1000)
    scv.tl.velocity(adata1, mode="dynamical")
    scv.tl.velocity_graph(adata1, n_jobs=n_jobs)
    scv.tl.latent_time(adata1)

elif mode == "deterministic":
    scv.tl.velocity(adata1, mode="deterministic")
    scv.tl.velocity_graph(adata1, n_jobs=n_jobs)
    # latent_time は dynamical 専用なので呼ばない

elif mode == "stochastic":
    scv.tl.velocity(adata1, mode="stochastic")
    scv.tl.velocity_graph(adata1, n_jobs=n_jobs)

else:
    raise SystemExit("mode must be 'stochastic' or 'deterministic' or 'dynamical'")

adata1.write(outfile, compression="lzf")
