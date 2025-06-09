# -*- coding: utf-8 -*-
import sys
import anndata
from scegot import scEGOT
import itertools
import plotly.io as pio
import pandas as pd
import numpy as np
import pickle

args = sys.argv
infile = args[1]
outfile = args[2]

# Set consstants
RANDOM_STATE = 2023
PCA_N_COMPONENTS = 30
GMM_CLUSTER_NUMBERS = [3, 3, 3, 3, 3]
UMAP_N_NEIGHBORS = 1000
DAY_NAMES = ["24h", "36h", "48h", "72h", "96h"]

# Loading
input_data = anndata.read_h5ad(infile)

# Create scEGOT instance
scegot = scEGOT(
    input_data,
    verbose=True,
    adata_day_key="cluster_day",
)

# Preprocessing
X, pca_model = scegot.preprocess(
    PCA_N_COMPONENTS,
    recode_params={},
    umi_target_sum=1e5,
    pca_random_state=RANDOM_STATE,
    pca_other_params={},
    apply_recode=True,
    apply_normalization_log1p=True,
    apply_normalization_umi=True,
    select_genes=True,
    n_select_genes=2000,
)

# Apply GMM
gmm_models, gmm_labels = scegot.fit_predict_gmm(
    n_components_list=GMM_CLUSTER_NUMBERS,
    covariance_type="full",
    max_iter=2000,
    n_init=10,
    random_state=RANDOM_STATE,
    gmm_other_params={},
)

# Apply UMAP
X_umap, umap_model = scegot.apply_umap(
    UMAP_N_NEIGHBORS,
    n_components=2,
    random_state=RANDOM_STATE,
    min_dist=0.8,
    umap_other_params={},
)

# Cell State Graph
cluster_names = scegot.generate_cluster_names_with_day()
G_pca = scegot.make_cell_state_graph(
    cluster_names,
    mode="pca",
    threshold=0.2,
)

G_umap = scegot.make_cell_state_graph(
    cluster_names,
    mode="umap",
    threshold=0.2,
)

# Velocity
velocities = scegot.calculate_cell_velocities()

# Waddington Potential
Wpotential, F_all = scegot.calculate_waddington_potential(
    n_neighbors=100,
    knn_other_params={},
)

# クラスタ番号と細胞型をDataFrameとしてまとめる
df = pd.DataFrame({
    "cluster": np.concatenate(scegot.gmm_labels),
    "celltype": input_data.obs["celltype"].values
})

# Output
bundle = {
    "input_data": input_data,
    "scegot": scegot,
    "X": X,
    "pca_model": pca_model,
    "gmm_models": gmm_models,
    "gmm_labels": gmm_labels,
    "X_umap": X_umap,
    "umap_model": umap_model,
    "cluster_names": cluster_names,
    "G_pca": G_pca,
    "G_umap": G_umap,
    "velocities": velocities,
    "Wpotential": Wpotential,
    "F_all": F_all,
    "df": df
}

# 書き出し（バイナリモード）
with open(outfile, "wb") as f:
    pickle.dump(bundle, f)
