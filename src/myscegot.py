# -*- coding: utf-8 -*-
import sys
import anndata
from scegot import scEGOT
import itertools
import plotly.io as pio

args = sys.argv
infile = args[1]
outfile = args[2]

# Set consstants
RANDOM_STATE = 2023
PCA_N_COMPONENTS = 30
GMM_CLUSTER_NUMBERS = [10, 10, 10, 10, 10]
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
G = scegot.make_cell_state_graph(
    cluster_names,
    mode="pca",
    threshold=0.2,
)

G_umap = scegot.make_cell_state_graph(
    cluster_names,
    mode="umap",
    threshold=0.2,
)

# Plot
scegot.plot_gmm_predictions(
    mode="pca",
    plot_gmm_means=True,
    figure_titles_without_gmm=[f"{name} on PCA" for name in DAY_NAMES],
    figure_titles_with_gmm=[f"{name} with GMM" for name in DAY_NAMES],
    cmap="plasma",
    save=True,
    save_paths=[f"./gmm_preds_pca_{day_name}.png" for day_name in DAY_NAMES]
)

scegot.plot_gmm_predictions(
    mode="umap",
    plot_gmm_means=True,
    figure_titles_without_gmm=[f"{name} on UMAP" for name in DAY_NAMES],
    figure_titles_with_gmm=[f"{name} with GMM" for name in DAY_NAMES],
    cmap="plasma",
    save=True,
    save_paths=[f"./gmm_preds_umap_{day_name}.png" for day_name in DAY_NAMES]
)

scegot.animatie_interpolated_distribution(
    cmap="gnuplot2",
    interpolate_interval=11,
    save=True,
    save_path="./cell_state_video.gif",
)

scegot.plot_simple_cell_state_graph(
    G_umap,
    layout="normal",
    save=True,
    save_path="./simple_cell_state_graph_umap.png"
)

scegot.plot_simple_cell_state_graph(
    G,
    layout="hierarchy",
    order="weight",
    save=True,
    save_path="./simple_cell_state_graph_hierarchy.png"
)

scegot.plot_true_and_interpolation_distributions(
    interpolate_index=2,
    mode="pca",
    n_samples=1000,
    t=0.5,
    plot_source_and_target=True,
    alpha_true=0.5,
    save=True,
    save_path="./true_and_interpolation_distributions_pca.png"
)

scegot.plot_true_and_interpolation_distributions(
    interpolate_index=2,
    mode="umap",
    n_samples=1000,
    t=0.5,
    plot_source_and_target=True,
    alpha_true=0.5,
    save=True,
    save_path="./true_and_interpolation_distributions_umap.png"
)

# Velocity
velocities = scegot.calculate_cell_velocities()
scegot.plot_cell_velocity(
    velocities,
    mode="pca",
    color_points="gmm",
    cluster_names=list(itertools.chain.from_iterable(cluster_names)),
    cmap="tab20",
    save=True,
    save_path="./cell_velocity_pca.png"
)

scegot.plot_cell_velocity(
    velocities,
    mode="umap",
    color_points="day",
    cmap="plasma",
    save=True,
    save_path="./cell_velocity_umap.png"
)

# Waddington Potential
Wpotential, F_all = scegot.calculate_waddington_potential(
    n_neighbors=100,
    knn_other_params={},
)

scegot.plot_waddington_potential(
    Wpotential,
    mode="pca",
    gene_name=None,
    save=True,
    save_path="./waddington_potential_pca_potential.html"
)

scegot.plot_waddington_potential(
    Wpotential,
    mode="umap",
    gene_name=None,
    save=True,
    save_path="./waddington_potential_umap_potential.html"
)

scegot.plot_waddington_potential_surface(
    Wpotential,
    mode="pca",
    save=True,
    save_path="./waddington_potential_surface_pca.html"
)

scegot.plot_waddington_potential_surface(
    Wpotential,
    mode="umap",
    save=True,
    save_path="./waddington_potential_surface_umap.html"
)