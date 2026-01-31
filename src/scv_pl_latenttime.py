# -*- coding: utf-8 -*-
import sys
import scvelo as scv
import scanpy
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

args = sys.argv
infile = args[1]
outfile = args[2]

# matplotlibの設定を最適化
matplotlib.use('Agg')  # 安定したバックエンド
plt.rcParams['savefig.dpi'] = 500
plt.rcParams['savefig.bbox'] = 'tight'
plt.rcParams['axes.facecolor'] = 'white'
plt.rcParams['savefig.facecolor'] = 'white'

# Parameters
scv.settings.verbosity = 3
scv.settings.presenter_view = True

# Loading
adata = scv.read(infile)

# プロット作成
fig, ax = plt.subplots(figsize=(10, 10), dpi=500)

# scatter plot with inferno colormap
scv.pl.scatter(adata, 
               color="latent_time", 
               basis='umap', 
               size=250,
               cmap='inferno',  # Infernoカラーマップ
               vmin=0,          # 最小値を0に固定
               vmax=1,          # 最大値を1に固定
               ax=ax,
               show=False,
               frameon=False)  # 枠を削除

# 軸とグリッドを削除してクリーンに
ax.set_xlabel('')
ax.set_ylabel('')
ax.set_title('Latent Time', fontsize=16, pad=20)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])

# スパイン（枠線）を削除
for spine in ax.spines.values():
    spine.set_visible(False)

# カラーバーの調整
cbar = ax.collections[0].colorbar
if cbar:
    cbar.ax.tick_params(labelsize=10)
    cbar.set_label('Latent Time', fontsize=12)
    # カラーバーの範囲も0-1に設定
    cbar.mappable.set_clim(vmin=0, vmax=1)

# 保存
plt.savefig(outfile, 
            dpi=500, 
            bbox_inches='tight', 
            facecolor='white',
            edgecolor='none')
plt.close()