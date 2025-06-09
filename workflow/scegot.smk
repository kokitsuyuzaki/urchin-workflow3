import pandas as pd
from snakemake.utils import min_version

GERM_LAYER = ['ectoderm', 'mesoderm', 'endoderm']

PLOT_FILES = ['gmm_preds_pca_24h.png', 'gmm_preds_pca_36h.png', 
              'gmm_preds_pca_48h.png', 'gmm_preds_pca_72h.png',
              'gmm_preds_pca_96h.png', 'gmm_preds_umap_24h.png',
              'gmm_preds_umap_36h.png', 'gmm_preds_umap_48h.png',
              'gmm_preds_umap_72h.png', 'gmm_preds_umap_96h.png',
              'cell_state_video.gif', 'simple_cell_state_graph_umap.png',
              'simple_cell_state_graph_hierarchy.png',
              'true_and_interpolation_distributions_pca.png',
              'true_and_interpolation_distributions_umap.png',
              'cell_velocity_pca_gmm.png',
              'cell_velocity_pca_day.png',
              'cell_velocity_umap_gmm.png',
              'cell_velocity_umap_day.png',
              'waddington_potential_pca_potential.html',
              'waddington_potential_umap_potential.html',
              'waddington_potential_surface_pca.html',
              'waddington_potential_surface_umap.html',
              'cluster_0_24h_celltype_pie.png',
              'cluster_0_36h_celltype_pie.png',
              'cluster_0_48h_celltype_pie.png',
              'cluster_0_72h_celltype_pie.png',
              'cluster_0_96h_celltype_pie.png',
              'cluster_1_24h_celltype_pie.png',
              'cluster_1_36h_celltype_pie.png',
              'cluster_1_48h_celltype_pie.png',
              'cluster_1_72h_celltype_pie.png',
              'cluster_1_96h_celltype_pie.png',
              'cluster_2_24h_celltype_pie.png',
              'cluster_2_36h_celltype_pie.png',
              'cluster_2_48h_celltype_pie.png',
              'cluster_2_72h_celltype_pie.png',
              'cluster_2_96h_celltype_pie.png']

#################################
# Setting
#################################
min_version("6.5.3")

rule all:
    input:
        expand('output/hpbase/cont/scegot_{germlayer}.pkl',
            germlayer=GERM_LAYER),
        expand('output/hpbase/DAPT/scegot_{germlayer}.pkl',
            germlayer=GERM_LAYER),
        expand('plot/hpbase/cont/scegot_{germlayer}/{file}',
            germlayer=GERM_LAYER, file=PLOT_FILES),
        expand('plot/hpbase/DAPT/scegot_{germlayer}/{file}',
            germlayer=GERM_LAYER, file=PLOT_FILES),

#################################
# Seurat => AnnData
#################################
rule seurat2anndata_cont_scegot:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/cont/seurat_scegot_ectoderm.h5ad',
        'output/hpbase/cont/seurat_scegot_mesoderm.h5ad',
        'output/hpbase/cont/seurat_scegot_endoderm.h5ad'
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
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/DAPT/seurat_scegot_ectoderm.h5ad',
        'output/hpbase/DAPT/seurat_scegot_mesoderm.h5ad',
        'output/hpbase/DAPT/seurat_scegot_endoderm.h5ad'
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
        'output/hpbase/cont/seurat_scegot_{germlayer}.h5ad'
    output:
        'output/hpbase/cont/scegot_{germlayer}.pkl'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scegot_cont_{germlayer}.txt'
    log:
        'logs/scegot_cont_{germlayer}.log'
    shell:
        'src/scegot.sh {input} {output} >& {log}'

rule scegot_DAPT:
    input:
        'output/hpbase/DAPT/seurat_scegot_{germlayer}.h5ad'
    output:
        'output/hpbase/DAPT/scegot_{germlayer}.pkl'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scegot_DAPT_{germlayer}.txt'
    log:
        'logs/scegot_DAPT_{germlayer}.log'
    shell:
        'src/scegot.sh {input} {output} >& {log}'


def aggregate_scegot_cont(germlayer):
    return expand('plot/hpbase/cont/scegot_{germlayer}/' + file,
                  germlayer=germlayer, file=PLOT_FILES)

def aggregate_scegot_DAPT(germlayer):
    return expand('plot/hpbase/DAPT/scegot_{germlayer}/' + file,
                  germlayer=germlayer, file=PLOT_FILES)

rule plot_scegot_cont:
    input:
        'output/hpbase/cont/scegot_{germlayer}.pkl'
    output:
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_pca_24h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_pca_36h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_pca_48h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_pca_72h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_pca_96h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_umap_24h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_umap_36h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_umap_48h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_umap_72h.png',
        'plot/hpbase/cont/scegot_{germlayer}/gmm_preds_umap_96h.png',
        'plot/hpbase/cont/scegot_{germlayer}/cell_state_video.gif',
        'plot/hpbase/cont/scegot_{germlayer}/simple_cell_state_graph_umap.png',
        'plot/hpbase/cont/scegot_{germlayer}/simple_cell_state_graph_hierarchy.png',
        'plot/hpbase/cont/scegot_{germlayer}/true_and_interpolation_distributions_pca.png',
        'plot/hpbase/cont/scegot_{germlayer}/true_and_interpolation_distributions_umap.png',
        'plot/hpbase/cont/scegot_{germlayer}/cell_velocity_pca_gmm.png',
        'plot/hpbase/cont/scegot_{germlayer}/cell_velocity_pca_day.png',
        'plot/hpbase/cont/scegot_{germlayer}/cell_velocity_umap_gmm.png',
        'plot/hpbase/cont/scegot_{germlayer}/cell_velocity_umap_day.png',
        'plot/hpbase/cont/scegot_{germlayer}/waddington_potential_pca_potential.html',
        'plot/hpbase/cont/scegot_{germlayer}/waddington_potential_umap_potential.html',
        'plot/hpbase/cont/scegot_{germlayer}/waddington_potential_surface_pca.html',
        'plot/hpbase/cont/scegot_{germlayer}/waddington_potential_surface_umap.html',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_0_24h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_0_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_0_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_0_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_0_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_1_24h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_1_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_1_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_1_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_1_96h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_2_24h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_2_36h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_2_48h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_2_72h_celltype_pie.png',
        'plot/hpbase/cont/scegot_{germlayer}/cluster_2_96h_celltype_pie.png'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_scegot_cont_{germlayer}.txt'
    log:
        'logs/plot_scegot_cont_{germlayer}.log'
    shell:
        'src/plot_scegot.sh {input} {output} >& {log}'

rule plot_scegot_DAPT:
    input:
        'output/hpbase/DAPT/scegot_{germlayer}.pkl'
    output:
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_pca_24h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_pca_36h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_pca_48h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_pca_72h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_pca_96h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_umap_24h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_umap_36h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_umap_48h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_umap_72h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/gmm_preds_umap_96h.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cell_state_video.gif',
        'plot/hpbase/DAPT/scegot_{germlayer}/simple_cell_state_graph_umap.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/simple_cell_state_graph_hierarchy.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/true_and_interpolation_distributions_pca.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/true_and_interpolation_distributions_umap.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cell_velocity_pca_gmm.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cell_velocity_pca_day.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cell_velocity_umap_gmm.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cell_velocity_umap_day.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/waddington_potential_pca_potential.html',
        'plot/hpbase/DAPT/scegot_{germlayer}/waddington_potential_umap_potential.html',
        'plot/hpbase/DAPT/scegot_{germlayer}/waddington_potential_surface_pca.html',
        'plot/hpbase/DAPT/scegot_{germlayer}/waddington_potential_surface_umap.html',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_0_24h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_0_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_0_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_0_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_0_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_1_24h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_1_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_1_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_1_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_1_96h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_2_24h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_2_36h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_2_48h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_2_72h_celltype_pie.png',
        'plot/hpbase/DAPT/scegot_{germlayer}/cluster_2_96h_celltype_pie.png'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_scegot_DAPT_{germlayer}.txt'
    log:
        'logs/plot_scegot_DAPT_{germlayer}.log'
    shell:
        'src/plot_scegot.sh {input} {output} >& {log}'
