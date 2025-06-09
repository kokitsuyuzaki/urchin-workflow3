# -*- coding: utf-8 -*-
import sys
import anndata
from scegot import scEGOT
import itertools
import plotly.io as pio
import matplotlib.pyplot as plt
from matplotlib import cm
import pandas as pd
import numpy as np
import pickle

args = sys.argv
infile = args[1]
outfiles = args[2:40]

# Set consstants
RANDOM_STATE = 2023
PCA_N_COMPONENTS = 30
GMM_CLUSTER_NUMBERS = [3, 3, 3, 3, 3]
UMAP_N_NEIGHBORS = 1000
DAY_NAMES = ["24h", "36h", "48h", "72h", "96h"]

# Loading
with open(infile, "rb") as f:
    loaded = pickle.load(f)

input_data = loaded["input_data"]
scegot = loaded["scegot"]
X = loaded["X"]
pca_model = loaded["pca_model"]
gmm_models = loaded["gmm_models"]
gmm_labels = loaded["gmm_labels"]
X_umap = loaded["X_umap"]
umap_model = loaded["umap_model"]
cluster_names = loaded["cluster_names"]
G_pca = loaded["G_pca"]
G_umap = loaded["G_umap"]
velocities = loaded["velocities"]
Wpotential = loaded["Wpotential"]
F_all = loaded["F_all"]
df = loaded["df"]

# Plot
scegot.plot_gmm_predictions(
    mode="pca",
    plot_gmm_means=True,
    figure_titles_without_gmm=[f"{name} on PCA" for name in DAY_NAMES],
    figure_titles_with_gmm=[f"{name} with GMM" for name in DAY_NAMES],
    cmap="plasma",
    save=True,
    save_paths=outfiles[:5]
)

scegot.plot_gmm_predictions(
    mode="umap",
    plot_gmm_means=True,
    figure_titles_without_gmm=[f"{name} on UMAP" for name in DAY_NAMES],
    figure_titles_with_gmm=[f"{name} with GMM" for name in DAY_NAMES],
    cmap="plasma",
    save=True,
    save_paths=outfiles[5:10]
)

scegot.animatie_interpolated_distribution(
    cmap="gnuplot2",
    interpolate_interval=11,
    save=True,
    save_path=outfiles[10],
)

scegot.plot_simple_cell_state_graph(
    G_umap,
    layout="normal",
    save=True,
    save_path=outfiles[11]
)

scegot.plot_simple_cell_state_graph(
    G_pca,
    layout="hierarchy",
    order="weight",
    save=True,
    save_path=outfiles[12]
)

scegot.plot_true_and_interpolation_distributions(
    interpolate_index=2,
    mode="pca",
    n_samples=1000,
    t=0.5,
    plot_source_and_target=True,
    alpha_true=0.5,
    save=True,
    save_path=outfiles[13]
)

scegot.plot_true_and_interpolation_distributions(
    interpolate_index=2,
    mode="umap",
    n_samples=1000,
    t=0.5,
    plot_source_and_target=True,
    alpha_true=0.5,
    save=True,
    save_path=outfiles[14]
)

# Velocity
scegot.plot_cell_velocity(
    velocities,
    mode="pca",
    color_points="gmm",
    cluster_names=list(itertools.chain.from_iterable(cluster_names)),
    cmap="tab20",
    save=True,
    save_path=outfiles[15]
)

scegot.plot_cell_velocity(
    velocities,
    mode="pca",
    color_points="day",
    cmap="tab20",
    save=True,
    save_path=outfiles[16]
)

scegot.plot_cell_velocity(
    velocities,
    mode="umap",
    color_points="gmm",
    cluster_names=list(itertools.chain.from_iterable(cluster_names)),
    cmap="tab20",
    save=True,
    save_path=outfiles[17]
)

scegot.plot_cell_velocity(
    velocities,
    mode="umap",
    color_points="day",
    cmap="tab20",
    save=True,
    save_path=outfiles[18]
)

# Waddington Potential
scegot.plot_waddington_potential(
    Wpotential,
    mode="pca",
    gene_name=None,
    save=True,
    save_path=outfiles[19]
)

scegot.plot_waddington_potential(
    Wpotential,
    mode="umap",
    gene_name=None,
    save=True,
    save_path=outfiles[20]
)

scegot.plot_waddington_potential_surface(
    Wpotential,
    mode="pca",
    save=True,
    save_path=outfiles[21]
)

scegot.plot_waddington_potential_surface(
    Wpotential,
    mode="umap",
    save=True,
    save_path=outfiles[22]
)

# 細胞型の固定順と色マップ
celltype_order = [
    "Aboral_ectoderm", "uncharacterized", "Ciliary_band",
    "Oral_ectoderm", "Neurons", "Endoderm",
    "Blastocoelar_cell", "Pigment", "Stomach",
    "Skeleton", "Stomach_Intestine", "Germ_line_future",
    "Intestine", "Pancreas", "Non_skeleton_mesoderm",
    "Anus", "Coelomic_pouch"
]
colors = plt.cm.tab20.colors
color_dict = {celltype: colors[i % len(colors)] for i, celltype in enumerate(celltype_order)}

# ファイル出力先
output_paths = outfiles[23:23 + 15]  # 3 clusters × 5 days = 15枚

# 時間軸ごとのオフセット
offset = 0
celltype_values = input_data.obs["celltype"].astype(str).values
n_clusters = 3  # 例：0, 1, 2 の3クラスタ

# すべての gmm_labels を day 単位で展開
gmm_labels_all = scegot.gmm_labels
day_names = scegot.day_names

# クラスタ → day という順でループ
out_index = 0
for cluster_id in range(n_clusters):
    offset = 0
    for day_idx, (labels_day, day_name) in enumerate(zip(gmm_labels_all, day_names)):
        n_cells = len(labels_day)
        celltypes_day = celltype_values[offset:offset + n_cells]
        labels_day = labels_day  # shape: (n_cells,)

        df = pd.DataFrame({
            "cluster": labels_day,
            "celltype": celltypes_day
        })

        subset = df[df["cluster"] == cluster_id]
        if subset.empty:
            out_index += 1
            continue

        counts = subset["celltype"].value_counts()
        labels = counts.index.tolist()
        pie_colors = [color_dict.get(label, "grey") for label in labels]

        plt.figure(figsize=(5, 5))
        plt.pie(counts, labels=labels, colors=pie_colors, autopct="%1.1f%%", startangle=90)
        plt.title(f"Cluster {cluster_id} - {day_name}")
        plt.tight_layout()
        plt.savefig(output_paths[out_index], transparent=True)
        plt.close()

        out_index += 1
        offset += n_cells