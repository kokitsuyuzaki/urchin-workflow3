import pandas as pd
from snakemake.utils import min_version

PLOT_FILES = [
    'gmm_preds_pca_36h.png',
    'gmm_preds_pca_48h.png',
    'gmm_preds_pca_72h.png',
    'gmm_preds_pca_96h.png',
    'gmm_preds_umap_36h.png',
    'gmm_preds_umap_48h.png',
    'gmm_preds_umap_72h.png',
    'gmm_preds_umap_96h.png',
    'cell_state_video.gif',
    'simple_cell_state_graph_umap.png',
    'simple_cell_state_graph_hierarchy.png',
    'true_and_interpolation_distributions_pca.png',
    'true_and_interpolation_distributions_umap.png',

    # --- Velocity ---
    'cell_velocity_pca_celltype.png',
    'cell_velocity_umap_celltype.png',
    'cell_velocity_pca_day.png',
    'cell_velocity_umap_day.png',
    'cell_velocity_pca_fieldonly.png',
    'cell_velocity_umap_fieldonly.png',

    'waddington_potential_pca_potential.html',
    'waddington_potential_umap_potential.html',
    'waddington_potential_surface_pca.html',
    'waddington_potential_surface_umap.html',

    'cluster_0_36h_celltype_pie.png',
    'cluster_0_48h_celltype_pie.png',
    'cluster_0_72h_celltype_pie.png',
    'cluster_0_96h_celltype_pie.png',
    'cluster_1_36h_celltype_pie.png',
    'cluster_1_48h_celltype_pie.png',
    'cluster_1_72h_celltype_pie.png',
    'cluster_1_96h_celltype_pie.png',
    'cluster_2_36h_celltype_pie.png',
    'cluster_2_48h_celltype_pie.png',
    'cluster_2_72h_celltype_pie.png',
    'cluster_2_96h_celltype_pie.png',
    'cluster_3_36h_celltype_pie.png',
    'cluster_3_48h_celltype_pie.png',
    'cluster_3_72h_celltype_pie.png',
    'cluster_3_96h_celltype_pie.png',
    'cluster_4_36h_celltype_pie.png',
    'cluster_4_48h_celltype_pie.png',
    'cluster_4_72h_celltype_pie.png',
    'cluster_4_96h_celltype_pie.png',
    'cluster_5_36h_celltype_pie.png',
    'cluster_5_48h_celltype_pie.png',
    'cluster_5_72h_celltype_pie.png',
    'cluster_5_96h_celltype_pie.png',
]

#################################
# Setting
#################################
min_version("6.5.3")

rule all:
    input:
        expand('plot/hpbase/cont/scegot_ectoderm/{file}',
            file=PLOT_FILES),
        expand('plot/hpbase/DAPT/scegot_ectoderm/{file}',
            file=PLOT_FILES),

#################################
# Seurat => AnnData
#################################
rule seurat2anndata_cont_scegot:
    input:
        'output/hpbase/integrated/seurat_annotated_landscaper.RData'
    output:
        'output/hpbase/cont/seurat_scegot_ectoderm.h5ad'
    container:
        'docker://koki/velocytor:20221015'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/seurat2anndata_cont_scegot.txt'
    log:
        'logs/seurat2anndata_cont_scegot.log'
    shell:
        'src/seurat2anndata_cont_scegot.sh {input} {output} >& {log}'

rule seurat2anndata_DAPT_scegot:
    input:
        'output/hpbase/integrated/seurat_annotated_landscaper.RData'
    output:
        'output/hpbase/DAPT/seurat_scegot_ectoderm.h5ad'
    container:
        'docker://koki/velocytor:20221015'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/seurat2anndata_DAPT_scegot.txt'
    log:
        'logs/seurat2anndata_DAPT_scegot.log'
    shell:
        'src/seurat2anndata_DAPT_scegot.sh {input} {output} >& {log}'

#################################
# scEGOT
#################################
rule scegot_cont:
    input:
        'output/hpbase/cont/seurat_scegot_ectoderm.h5ad'
    output:
        'output/hpbase/cont/scegot_ectoderm.pkl'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scegot_cont_ectoderm.txt'
    log:
        'logs/scegot_cont_ectoderm.log'
    shell:
        'src/scegot.sh {input} {output} >& {log}'

rule scegot_DAPT:
    input:
        'output/hpbase/DAPT/seurat_scegot_ectoderm.h5ad'
    output:
        'output/hpbase/DAPT/scegot_ectoderm.pkl'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scegot_DAPT_ectoderm.txt'
    log:
        'logs/scegot_DAPT_ectoderm.log'
    shell:
        'src/scegot.sh {input} {output} >& {log}'

rule plot_scegot_cont:
    input:
        'output/hpbase/cont/scegot_ectoderm.pkl'
    output:
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_pca_36h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_pca_48h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_pca_72h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_pca_96h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_umap_36h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_umap_48h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_umap_72h.png',
        'plot/hpbase/cont/scegot_ectoderm/gmm_preds_umap_96h.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_state_video.gif',
        'plot/hpbase/cont/scegot_ectoderm/simple_cell_state_graph_umap.png',
        'plot/hpbase/cont/scegot_ectoderm/simple_cell_state_graph_hierarchy.png',
        'plot/hpbase/cont/scegot_ectoderm/true_and_interpolation_distributions_pca.png',
        'plot/hpbase/cont/scegot_ectoderm/true_and_interpolation_distributions_umap.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_velocity_pca_celltype.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_velocity_umap_celltype.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_velocity_pca_day.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_velocity_umap_day.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_velocity_pca_fieldonly.png',
        'plot/hpbase/cont/scegot_ectoderm/cell_velocity_umap_fieldonly.png',
        'plot/hpbase/cont/scegot_ectoderm/waddington_potential_pca_potential.html',
        'plot/hpbase/cont/scegot_ectoderm/waddington_potential_umap_potential.html',
        'plot/hpbase/cont/scegot_ectoderm/waddington_potential_surface_pca.html',
        'plot/hpbase/cont/scegot_ectoderm/waddington_potential_surface_umap.html',
        'plot/hpbase/cont/scegot_ectoderm/cluster_0_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_0_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_0_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_0_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_1_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_1_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_1_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_1_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_2_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_2_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_2_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_2_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_3_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_3_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_3_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_3_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_4_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_4_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_4_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_4_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_5_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_5_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_5_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_ectoderm/cluster_5_96h_celltype_pie.png'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_scegot_cont_ectoderm.txt'
    log:
        'logs/plot_scegot_cont_ectoderm.log'
    shell:
        # src/plot_scegot.sh は「python plot_scegot.py {input} {output}...」形式で
        # outputをそのまま全部渡している前提
        'src/plot_scegot.sh {input} {output} >& {log}'


rule plot_scegot_DAPT:
    input:
        'output/hpbase/DAPT/scegot_ectoderm.pkl'
    output:
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_pca_36h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_pca_48h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_pca_72h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_pca_96h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_umap_36h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_umap_48h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_umap_72h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/gmm_preds_umap_96h.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_state_video.gif',
        'plot/hpbase/DAPT/scegot_ectoderm/simple_cell_state_graph_umap.png',
        'plot/hpbase/DAPT/scegot_ectoderm/simple_cell_state_graph_hierarchy.png',
        'plot/hpbase/DAPT/scegot_ectoderm/true_and_interpolation_distributions_pca.png',
        'plot/hpbase/DAPT/scegot_ectoderm/true_and_interpolation_distributions_umap.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_velocity_pca_celltype.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_velocity_umap_celltype.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_velocity_pca_day.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_velocity_umap_day.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_velocity_pca_fieldonly.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cell_velocity_umap_fieldonly.png',
        'plot/hpbase/DAPT/scegot_ectoderm/waddington_potential_pca_potential.html',
        'plot/hpbase/DAPT/scegot_ectoderm/waddington_potential_umap_potential.html',
        'plot/hpbase/DAPT/scegot_ectoderm/waddington_potential_surface_pca.html',
        'plot/hpbase/DAPT/scegot_ectoderm/waddington_potential_surface_umap.html',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_0_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_0_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_0_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_0_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_1_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_1_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_1_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_1_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_2_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_2_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_2_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_2_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_3_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_3_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_3_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_3_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_4_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_4_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_4_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_4_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_5_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_5_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_5_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_ectoderm/cluster_5_96h_celltype_pie.png'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_scegot_DAPT_ectoderm.txt'
    log:
        'logs/plot_scegot_DAPT_ectoderm.log'
    shell:
        'src/plot_scegot.sh {input} {output} >& {log}'
