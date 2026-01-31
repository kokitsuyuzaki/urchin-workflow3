# -*- coding: utf-8 -*-
import sys
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PathCollection

infile = sys.argv[1]
outfile = sys.argv[2]

adata = sc.read_h5ad(infile)

# ========= celltype / celltype_colors =========
if "celltype" in adata.obs and "celltype_colors" in adata.obs:
    celltype_color_map = {}
    for ct in adata.obs["celltype"].unique():
        color = adata.obs.loc[adata.obs["celltype"] == ct, "celltype_colors"].iloc[0]
        celltype_color_map[ct] = color

    celltype_order = sorted(celltype_color_map.keys())
    adata.obs["celltype"] = pd.Categorical(
        adata.obs["celltype"], categories=celltype_order
    )
    adata.uns["celltype_colors"] = [celltype_color_map[ct] for ct in celltype_order]

# ========= neighbors / PAGA / UMAP =========
if "neighbors" not in adata.uns:
    sc.pp.neighbors(adata)

group_key = "celltype" if "celltype" in adata.obs else "seurat_clusters"
adata.obs[group_key] = adata.obs[group_key].astype("category")

if "paga" not in adata.uns:
    sc.tl.paga(adata, groups=group_key)

if "X_umap" not in adata.obsm:
    sc.tl.umap(adata)

# ========= PAGA node positions = UMAP重心 =========
umap = adata.obsm["X_umap"][:, :2]
cats = adata.obs[group_key].cat.categories

pos = np.zeros((len(cats), 2))
for i, ct in enumerate(cats):
    idx = (adata.obs[group_key] == ct).values
    pos[i] = umap[idx].mean(axis=0)

# ========= Plot =========
sc.set_figure_params(dpi=500, frameon=False)
fig, ax = plt.subplots(figsize=(10, 10))

# --- UMAP points ---
sc.pl.umap(
    adata,
    color=group_key,
    size=80,          # データ点は大きめ
    alpha=0.85,
    show=False,
    ax=ax,
)

# legend / title を消す
ax.set_title("")
leg = ax.get_legend()
if leg is not None:
    leg.remove()

# --- PAGA描画前の状態を保存 ---
before_cols = list(ax.collections)
before_texts = list(ax.texts)

# --- PAGA抽象グラフ ---
sc.pl.paga(
    adata,
    color=group_key,
    threshold=0.3,
    pos=pos,
    show=False,
    ax=ax,
)

# ========= 後処理 =========

# 1) エッジ：元の太さ ×3
for c in ax.collections:
    if isinstance(c, LineCollection) and c not in before_cols:
        lw = np.asarray(c.get_linewidths(), dtype=float)
        c.set_linewidth(lw * 3.0)
        c.set_color("#000000")

# 2) PAGAノード：黒縁を付ける
for c in ax.collections:
    if isinstance(c, PathCollection) and c not in before_cols:
        sz = np.asarray(c.get_sizes(), dtype=float)
        c.set_sizes(sz * 3.5)        # ← ここを好みで調整（例: 2.0, 3.0）
        c.set_edgecolor("black")
        c.set_linewidth(1.5)   # 縁の太さ（好みで調整）
        c.set_zorder(5)        # 前面に出す

# 3) クラスタ名ラベルを完全に消す
for t in ax.texts:
    if t not in before_texts:
        t.set_visible(False)

fig.savefig(outfile, dpi=500, bbox_inches="tight")
plt.close(fig)
