# -*- coding: utf-8 -*-
import sys
import anndata
from scegot import scEGOT
import itertools
import plotly.io as pio
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import pickle

# 追加
import seaborn as sns

args = sys.argv
infile = args[1]
outfiles = args[2:]

# Set constants
RANDOM_STATE = 2023
PCA_N_COMPONENTS = 30
GMM_CLUSTER_NUMBERS = [6, 6, 6]  # ここは生成側で使うのでこのまま
UMAP_N_NEIGHBORS = 1000
DAY_NAMES = ["36h", "48h", "72h", "96h"]  # 4日

# day色（固定）
DAY_COLOR = {
    "36h": "#3498DB",  # 青系
    "48h": "#27AE60",  # 緑系
    "72h": "#F39C12",  # オレンジ系
    "96h": "#8E44AD",  # 紫系
}

# celltype色（必要ならここを編集）
celltype_order = [
    "Aboral_ectoderm", "uncharacterized", "Ciliary_band",
    "Oral_ectoderm", "Neurons", "Anus",
]
CELLTYPE_COLOR = {
    "Aboral_ectoderm": "#008080",
    "Oral_ectoderm":   "#FFFF00",
    "Ciliary_band":    "#00FF26",
    "Neurons":         "#FF00FF",
    "Anus":            "#3F3F7F",
    "uncharacterized": "#7f7f7f",
}

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
velocities = loaded["velocities"]          # 多くの場合 (18335, 30)（96hは含まれない仕様）
Wpotential = loaded["Wpotential"]
F_all = loaded["F_all"]
df = loaded["df"]


# -----------------------------
# helpers (plots)
# -----------------------------
def _save_transparent(fig, path, dpi=300):
    fig.savefig(path, dpi=dpi, transparent=True)
    plt.close(fig)


def plot_vector_field_only(scegot, velocities, mode, save_path,
                           density=2, arrow_size=1, linewidth=1.5):
    """
    ベクトル場のみ（点なし）、背景透明。
    velocities の行数に合わせて X_{pca,umap} を切って作図する。
    """
    import scanpy as sc
    import scvelo as scv

    if mode not in ["pca", "umap"]:
        raise ValueError("mode must be 'pca' or 'umap'")

    # velocity がある分だけに合わせる（96hがないなら自動で落ちる）
    X_pca_all = pd.concat(scegot.X_pca)
    X_umap_all = pd.concat(scegot.X_umap)
    nV = velocities.shape[0]

    X_pca = X_pca_all.iloc[:nV, :]

    adata = anndata.AnnData(
        X_pca.values,
        obsm={"X_pca": X_pca.values},
        layers={"velocity": velocities.values, "spliced": X_pca.values},
    )
    if mode == "umap":
        X_umap = X_umap_all.iloc[:nV, :]
        adata.obsm["X_umap"] = X_umap.values

    # neighborsはPCAで構築（scEGOTのplot_cell_velocity踏襲）
    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=scegot.pca_model.n_components_)
    scv.tl.velocity_graph(adata)

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True, facecolor="none")
    ax.set_facecolor("none")

    # streamのみ（点なし）
    scv.pl.velocity_embedding_stream(
        adata,
        basis=mode,
        color="black",
        vkey="velocity",
        title="",
        density=density,
        alpha=0.0,           # 点なし
        fontsize=14,
        legend_fontsize=0,
        legend_loc=None,
        arrow_size=arrow_size,
        linewidth=linewidth,
        ax=ax,
        show=False,
        X_grid=None,
        V_grid=None,
        sort_order=True,
        size=50,
        colorbar=False,
    )

    ax.axis("off")
    _save_transparent(fig, save_path)


def plot_points_only_celltype(scegot, input_data, mode, save_path,
                              point_alpha=0.8, size_points=10):
    """
    cell_velocity_{pca,umap}_celltype.png:
      点のみ（ベクトル場なし）、celltypeで色付け。
      ※velocitiesに依存しない（=96hも描ける）
    """
    if mode not in ["pca", "umap"]:
        raise ValueError("mode must be 'pca' or 'umap'")

    X_list = scegot.X_pca if mode == "pca" else scegot.X_umap
    X_cat = pd.concat(X_list).iloc[:, :2]

    celltype = input_data.obs["celltype"].astype(str).str.strip().values
    if len(celltype) != len(X_cat):
        raise ValueError(f"Length mismatch: celltype={len(celltype)} vs X_{mode}={len(X_cat)}")

    colors = [CELLTYPE_COLOR.get(ct, "#999999") for ct in celltype]

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True, facecolor="none")
    ax.set_facecolor("none")

    ax.scatter(
        X_cat.iloc[:, 0],
        X_cat.iloc[:, 1],
        c=colors,
        s=size_points,
        alpha=point_alpha,
        edgecolors="none",
    )

    ax.axis("off")
    _save_transparent(fig, save_path)


def plot_points_only_day(scegot, mode, save_path,
                         point_alpha=0.8, size_points=10):
    """
    cell_velocity_{pca,umap}_day.png:
      点のみ（ベクトル場なし）、dayで色付け。
      ※velocitiesに依存しない（=96hも描ける）
    """
    if mode not in ["pca", "umap"]:
        raise ValueError("mode must be 'pca' or 'umap'")

    day_names = list(scegot.day_names)
    X_list = scegot.X_pca if mode == "pca" else scegot.X_umap
    X_cat = pd.concat(X_list).iloc[:, :2]

    # dayごとに同数の色を繰り返し
    color_vec = []
    for i, d in enumerate(day_names):
        color_vec += [DAY_COLOR.get(d, "#777777")] * len(X_list[i])

    if len(color_vec) != len(X_cat):
        raise ValueError(f"Length mismatch: colors={len(color_vec)} vs X_{mode}={len(X_cat)}")

    sns.set_style("white")
    fig, ax = plt.subplots(figsize=(8, 6), tight_layout=True, facecolor="none")
    ax.set_facecolor("none")

    ax.scatter(
        X_cat.iloc[:, 0],
        X_cat.iloc[:, 1],
        c=color_vec,
        s=size_points,
        alpha=point_alpha,
        edgecolors="none",
    )

    ax.axis("off")
    _save_transparent(fig, save_path)


# -----------------------------
# Main pipeline (same as before)
# -----------------------------

# ---- GMM predictions (PCA/UMAP) ----
scegot.plot_gmm_predictions(
    mode="pca",
    plot_gmm_means=True,
    figure_titles_without_gmm=[f"{name} on PCA" for name in DAY_NAMES],
    figure_titles_with_gmm=[f"{name} with GMM" for name in DAY_NAMES],
    cmap="plasma",
    save=True,
    save_paths=outfiles[0:4],   # gmm_preds_pca_36/48/72/96h.png
)

scegot.plot_gmm_predictions(
    mode="umap",
    plot_gmm_means=True,
    figure_titles_without_gmm=[f"{name} on UMAP" for name in DAY_NAMES],
    figure_titles_with_gmm=[f"{name} with GMM" for name in DAY_NAMES],
    cmap="plasma",
    save=True,
    save_paths=outfiles[4:8],   # gmm_preds_umap_36/48/72/96h.png
)

# ---- animated interpolation ----
scegot.animatie_interpolated_distribution(
    cmap="gnuplot2",
    interpolate_interval=11,
    save=True,
    save_path=outfiles[8],      # cell_state_video.gif
)

# ---- cell state graphs ----
scegot.plot_simple_cell_state_graph(
    G_umap,
    layout="normal",
    save=True,
    save_path=outfiles[9],      # simple_cell_state_graph_umap.png
)

scegot.plot_simple_cell_state_graph(
    G_pca,
    layout="hierarchy",
    order="weight",
    save=True,
    save_path=outfiles[10],     # simple_cell_state_graph_hierarchy.png
)

# ---- interpolation distributions ----
scegot.plot_true_and_interpolation_distributions(
    interpolate_index=0,
    mode="pca",
    n_samples=1000,
    t=0.5,
    plot_source_and_target=True,
    alpha_true=0.5,
    save=True,
    save_path=outfiles[11],     # true_and_interpolation_distributions_pca.png
)

scegot.plot_true_and_interpolation_distributions(
    interpolate_index=1,
    mode="umap",
    n_samples=1000,
    t=0.5,
    plot_source_and_target=True,
    alpha_true=0.5,
    save=True,
    save_path=outfiles[12],     # true_and_interpolation_distributions_umap.png
)

# ---- Velocity outputs (new spec) ----
# index mapping (Snakemakeで合わせた前提):
# 13: cell_velocity_pca_celltype.png
# 14: cell_velocity_umap_celltype.png
# 15: cell_velocity_pca_day.png
# 16: cell_velocity_umap_day.png
# 17: cell_velocity_pca_fieldonly.png
# 18: cell_velocity_umap_fieldonly.png

# 点のみ（celltype）: velocities不要（=96hも描ける）
plot_points_only_celltype(
    scegot=scegot,
    input_data=input_data,
    mode="pca",
    save_path=outfiles[13],
    point_alpha=0.8,
    size_points=10,
)
plot_points_only_celltype(
    scegot=scegot,
    input_data=input_data,
    mode="umap",
    save_path=outfiles[14],
    point_alpha=0.8,
    size_points=10,
)

# 点のみ（day）: velocities不要（=96hも描ける）
plot_points_only_day(
    scegot=scegot,
    mode="pca",
    save_path=outfiles[15],
    point_alpha=0.8,
    size_points=10,
)
plot_points_only_day(
    scegot=scegot,
    mode="umap",
    save_path=outfiles[16],
    point_alpha=0.8,
    size_points=10,
)

# ベクトル場のみ（点なし、背景透明）: velocitiesを使う（=96hは無い仕様のまま）
plot_vector_field_only(
    scegot=scegot,
    velocities=velocities,
    mode="pca",
    save_path=outfiles[17],
)
plot_vector_field_only(
    scegot=scegot,
    velocities=velocities,
    mode="umap",
    save_path=outfiles[18],
)

# ---- Waddington Potential ----
scegot.plot_waddington_potential(
    Wpotential,
    mode="pca",
    gene_name=None,
    save=True,
    save_path=outfiles[19],     # waddington_potential_pca_potential.html
)

scegot.plot_waddington_potential(
    Wpotential,
    mode="umap",
    gene_name=None,
    save=True,
    save_path=outfiles[20],     # waddington_potential_umap_potential.html
)

scegot.plot_waddington_potential_surface(
    Wpotential,
    mode="pca",
    save=True,
    save_path=outfiles[21],     # waddington_potential_surface_pca.html
)

scegot.plot_waddington_potential_surface(
    Wpotential,
    mode="umap",
    save=True,
    save_path=outfiles[22],     # waddington_potential_surface_umap.html
)

# ---- celltype pie charts (clusters × days) ----
colors = plt.cm.tab20.colors
color_dict = {celltype: colors[i % len(colors)] for i, celltype in enumerate(celltype_order)}

# 6 clusters × 4 days = 24枚 → outfiles[23:47]
output_paths = outfiles[23:23 + 24]
celltype_values = input_data.obs["celltype"].astype(str).values
n_clusters = 6

gmm_labels_all = scegot.gmm_labels
day_names = ["36h", "48h", "72h", "96h"]

# 各日のセル数と累積オフセット
day_sizes = [len(lbl) for lbl in gmm_labels_all]
cum = np.cumsum([0] + day_sizes)

out_index = 0
for cluster_id in range(n_clusters):       # 0..5
    for day_idx, (labels_day, day_name) in enumerate(zip(gmm_labels_all, day_names)):
        start, end = cum[day_idx], cum[day_idx + 1]
        celltypes_day = celltype_values[start:end]

        df_plot = pd.DataFrame({"cluster": labels_day, "celltype": celltypes_day})
        subset = df_plot[df_plot["cluster"] == cluster_id]

        plt.figure(figsize=(5, 5))
        if subset.empty:
            plt.title(f"Cluster {cluster_id} - {day_name} (empty)")
        else:
            counts = subset["celltype"].value_counts()
            labels = counts.index.tolist()
            pie_colors = [color_dict.get(label, "grey") for label in labels]
            plt.pie(
                counts,
                labels=labels,
                colors=pie_colors,
                autopct="%1.1f%%",
                startangle=90,
            )
            plt.title(f"Cluster {cluster_id} - {day_name}")
        plt.tight_layout()
        plt.savefig(output_paths[out_index], transparent=True)
        plt.close()

        out_index += 1
