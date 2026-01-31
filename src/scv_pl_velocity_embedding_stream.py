# -*- coding: utf-8 -*-
import sys
import scvelo as scv
import scanpy
import pandas as pd
import numpy as np

args = sys.argv
infile = args[1]
outfile = args[2]

# Parameters
scv.settings.verbosity = 3
scv.settings.presenter_view = True
scv.set_figure_params('scvelo')

# Loading
adata = scv.read(infile)

print("Available columns in adata.obs:")
print(adata.obs.columns.tolist())

# celltype_colorsからカラーマップを作成
if 'celltype' in adata.obs and 'celltype_colors' in adata.obs:
    # ユニークなcelltypeと対応する色を取得
    celltype_color_map = {}
    for celltype in adata.obs['celltype'].unique():
        # そのcelltypeの最初の色を取得
        color = adata.obs[adata.obs['celltype'] == celltype]['celltype_colors'].iloc[0]
        celltype_color_map[celltype] = color
    
    # 色のリストを作成（celltypeのカテゴリ順）
    celltype_order = sorted(celltype_color_map.keys())
    palette = [celltype_color_map[ct] for ct in celltype_order]
    
    print("Celltype color mapping:")
    for ct, color in celltype_color_map.items():
        print(f"  {ct}: {color}")
    
    # unsにも保存（scveloが参照する）
    adata.uns['celltype_colors'] = palette
    
    # Embedding with custom colors
    scv.pl.velocity_embedding_stream(
        adata, 
        basis='umap', 
        save=outfile, 
        dpi=500, 
        figsize=(10,10),
        color="celltype",
        palette=celltype_color_map,
        linewidth=2.5,
        size=80,
        alpha=0.8,
        density=0.5,        # 1.5 → 0.5 に減らす（流線を減らす）
        arrow_color='black',
        arrow_size=3.0,
        title='',           # タイトルを消す
        legend_loc='none'   # 凡例（細胞型ラベル）を消す
    )
else:
    # celltype_colorsがない場合はデフォルト
    scv.pl.velocity_embedding_stream(
        adata, 
        basis='umap', 
        save=outfile, 
        dpi=500, 
        figsize=(10,10),
        color="celltype",
        linewidth=2.5,
        size=80,
        density=0.5,        # 流線を減らす
        title='',           # タイトルを消す
        legend_loc='none'   # 凡例を消す
    )