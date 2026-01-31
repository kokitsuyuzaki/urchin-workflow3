import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

METHODS = ["paga", "paga_tree", "comp1", "mst", "angle"]
MODES = ['deterministic', 'stochastic', 'dynamical']
INTEGRATED_PLOTS = ['elbowplot.png', 'barplot.png',
    'dimplot_cluster.png', 'dimplot_cluster_splitby.png',
    'featureplot_ncount_rna.png', 'featureplot_ncount_rna_splitby.png',
    'featureplot_nfeature_rna.png', 'featureplot_nfeature_rna_splitby.png',
    'featureplot_percent_mt.png', 'featureplot_percent_mt_splitby.png',
    'featureplot_percent_rb.png', 'featureplot_percent_rb_splitby.png',
    'ridgeplot_cellcycle.png', 'dimplot_cellcycle.png',
    'dimplot_cellcycle_splitby.png',
    'marker/FINISH_marker', 'marker/FINISH_cluster_marker',
    'featureplot_doublet.png', 'featureplot_doublet_splitby.png',
    'plot_cells_trajectory.png']

CONT_PLOTS = ['elbowplot.png', 'dimplot_cluster.png', 'dimplot_cluster_splitby.png']

DAPT_PLOTS = CONT_PLOTS

LABELTYPES = ['clusters', 'sample', 'celltype']

rule all:
    input:
        # Basic Plots
        expand('plot/hpbase/cont_stratified/dynverse_{method}_celltype_{labeltype}.png',
            method=METHODS, labeltype=LABELTYPES),
        expand('plot/hpbase/DAPT_stratified/dynverse_{method}_celltype_{labeltype}.png',
            method=METHODS, labeltype=LABELTYPES),
        expand('plot/hpbase/cont/{cont_plot}',
            cont_plot=CONT_PLOTS),
        expand('plot/hpbase/DAPT/{dapt_plot}',
            dapt_plot=DAPT_PLOTS),
        'plot/hpbase/cont_stratified/plot_cells_trajectory.png',
        'plot/hpbase/DAPT_stratified/plot_cells_trajectory.png',
        # Control
        'plot/hpbase/cont/scv_pl_proportions_cont_dynamical.png',
        expand('plot/hpbase/cont/scv_pl_velocity_embedding_stream_cont_{mode}.png',
            mode=MODES),
        expand('plot/hpbase/cont/scv_pl_velocity_embedding_cont_{mode}.png',
            mode=MODES),
        'plot/hpbase/cont/scv_pl_latenttime_cont_dynamical.png',
        'plot/hpbase/cont/scv_pl_heatmap_cont_dynamical.png',
        'plot/hpbase/cont/scv_pl_velocity_markers_cont_dynamical_1.png',
        'plot/hpbase/cont/scv_pl_velocity_markers_cont_dynamical_2.png',
        'plot/hpbase/cont/scv_pl_velocity_markers_cont_dynamical_3.png',
        'plot/hpbase/cont/scv_pl_paga_cont.png',
        # DAPT
        'plot/hpbase/DAPT/scv_pl_proportions_DAPT_dynamical.png',
        expand('plot/hpbase/DAPT/scv_pl_velocity_embedding_stream_DAPT_{mode}.png',
            mode=MODES),
        expand('plot/hpbase/DAPT/scv_pl_velocity_embedding_DAPT_{mode}.png',
            mode=MODES),
        'plot/hpbase/DAPT/scv_pl_latenttime_DAPT_dynamical.png',
        'plot/hpbase/DAPT/scv_pl_heatmap_DAPT_dynamical.png',
        'plot/hpbase/DAPT/scv_pl_velocity_markers_DAPT_dynamical_1.png',
        'plot/hpbase/DAPT/scv_pl_velocity_markers_DAPT_dynamical_2.png',
        'plot/hpbase/DAPT/scv_pl_velocity_markers_DAPT_dynamical_3.png',
        'plot/hpbase/DAPT/scv_pl_paga_DAPT.png',

#################################
# Elbow Plot
#################################
rule elbowplot_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'plot/hpbase/integrated/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/elbowplot_integrated_integrated.txt'
    log:
        'logs/elbowplot_integrated_integrated.log'
    shell:
        'src/elbowplot_integrated.sh {input} {output} >& {log}'

rule elbowplot_cont:
    input:
        'output/hpbase/cont/seurat.RData'
    output:
        'plot/hpbase/cont/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/elbowplot_cont.txt'
    log:
        'logs/elbowplot_cont.log'
    shell:
        'src/elbowplot_integrated.sh {input} {output} >& {log}'

rule elbowplot_DAPT:
    input:
        'output/hpbase/DAPT/seurat.RData'
    output:
        'plot/hpbase/DAPT/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/elbowplot_DAPT.txt'
    log:
        'logs/elbowplot_DAPT.log'
    shell:
        'src/elbowplot_integrated.sh {input} {output} >& {log}'

#################################
# Barplot
#################################
rule barplot_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'plot/hpbase/integrated/barplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/barplot_integrated_integrated.txt'
    log:
        'logs/barplot_integrated_integrated.log'
    shell:
        'src/barplot_integrated.sh {input} {output} >& {log}'

#################################
# Cluster Label
#################################
rule dimplot_cluster_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'plot/hpbase/integrated/dimplot_cluster.png',
        'plot/hpbase/integrated/dimplot_cluster_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_cluster_integrated_integrated.txt'
    log:
        'logs/dimplot_cluster_integrated_integrated.log'
    shell:
        'src/dimplot_cluster_integrated.sh {input} {output} >& {log}'

rule dimplot_cluster_cont:
    input:
        'output/hpbase/cont/seurat.RData'
    output:
        'plot/hpbase/cont/dimplot_cluster.png',
        'plot/hpbase/cont/dimplot_cluster_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_cluster_cont.txt'
    log:
        'logs/dimplot_cluster_cont.log'
    shell:
        'src/dimplot_cluster_integrated.sh {input} {output} >& {log}'

rule dimplot_cluster_DAPT:
    input:
        'output/hpbase/DAPT/seurat.RData'
    output:
        'plot/hpbase/DAPT/dimplot_cluster.png',
        'plot/hpbase/DAPT/dimplot_cluster_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_cluster_DAPT.txt'
    log:
        'logs/dimplot_cluster_DAPT.log'
    shell:
        'src/dimplot_cluster_integrated.sh {input} {output} >& {log}'

#################################
# Number of RNA counts
#################################
rule featureplot_ncount_rna_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'plot/hpbase/integrated/featureplot_ncount_rna.png',
        'plot/hpbase/integrated/featureplot_ncount_rna_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_ncount_rna_integrated_integrated.txt'
    log:
        'logs/featureplot_ncount_rna_integrated_integrated.log'
    shell:
        'src/featureplot_ncount_rna_integrated.sh {input} {output} >& {log}'

#################################
# Number of detected RNAs
#################################
rule featureplot_nfeature_rna_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'plot/hpbase/integrated/featureplot_nfeature_rna.png',
        'plot/hpbase/integrated/featureplot_nfeature_rna_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_nfeature_rna_integrated_integrated.txt'
    log:
        'logs/featureplot_nfeature_rna_integrated_integrated.log'
    shell:
        'src/featureplot_nfeature_rna_integrated.sh {input} {output} >& {log}'

#################################
# Percentage of mitochondria genes' expression
#################################
rule featureplot_percent_mt_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/hpbase/integrated/featureplot_percent_mt.png',
        'plot/hpbase/integrated/featureplot_percent_mt_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_percent_mt_integrated_integrated.txt'
    log:
        'logs/featureplot_percent_mt_integrated_integrated.log'
    shell:
        'src/featureplot_percent_mt_integrated.sh {input} {output} >& {log}'

#################################
# Percentage of ribosome genes' expression
#################################
rule featureplot_percent_rb_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/hpbase/integrated/featureplot_percent_rb.png',
        'plot/hpbase/integrated/featureplot_percent_rb_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_percent_rb_integrated_integrated.txt'
    log:
        'logs/featureplot_percent_rb_integrated_integrated.log'
    shell:
        'src/featureplot_percent_rb_integrated.sh {input} {output} >& {log}'

#################################
# Cell cycle score
#################################
rule dimplot_cellcycle_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/hpbase/integrated/ridgeplot_cellcycle.png',
        'plot/hpbase/integrated/dimplot_cellcycle.png',
        'plot/hpbase/integrated/dimplot_cellcycle_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_cellcycle_integrated_integrated.txt'
    log:
        'logs/dimplot_cellcycle_integrated_integrated.log'
    shell:
        'src/dimplot_cellcycle_integrated.sh {input} {output} >& {log}'

#################################
# Marker Genes (Prepared)
#################################
rule featureplot_marker_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/marker.RData'
    output:
        'plot/hpbase/integrated/marker/FINISH_marker'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_marker_integrated_integrated.txt'
    log:
        'logs/featureplot_marker_integrated_integrated.log'
    shell:
        'src/featureplot_marker_integrated.sh {input} {output} >& {log}'

#################################
# Marker Genes (Detected in each Cluster)
#################################
rule featureplot_cluster_marker_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'output/hpbase/integrated/markers.xlsx'
    output:
        'plot/hpbase/integrated/marker/FINISH_cluster_marker'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_cluster_marker_integrated_integrated.txt'
    log:
        'logs/featureplot_cluster_marker_integrated_integrated.log'
    shell:
        'src/featureplot_cluster_marker_integrated.sh {input} {output} >& {log}'

#################################
# Doublet Density
#################################
rule featureplot_doublet_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'output/hpbase/integrated/scdblfinder.RData'
    output:
        'plot/hpbase/integrated/featureplot_doublet.png',
        'plot/hpbase/integrated/featureplot_doublet_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_doublet_integrated_integrated.txt'
    log:
        'logs/featureplot_doublet_integrated_integrated.log'
    shell:
        'src/featureplot_doublet_integrated.sh {input} {output} >& {log}'

#################################
# Trajectory Inference
#################################
rule plot_cells_trajectory_cont:
    input:
        'output/hpbase/cont_stratified/monocle3.RData'
    output:
        'plot/hpbase/cont_stratified/plot_cells_trajectory.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_cells_trajectory_cont.txt'
    log:
        'logs/plot_cells_trajectory_cont.log'
    shell:
        'src/plot_cells_trajectory.sh {input} {output} >& {log}'

rule plot_cells_trajectory_DAPT:
    input:
        'output/hpbase/DAPT_stratified/monocle3.RData'
    output:
        'plot/hpbase/DAPT_stratified/plot_cells_trajectory.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_cells_trajectory_DAPT.txt'
    log:
        'logs/plot_cells_trajectory_DAPT.log'
    shell:
        'src/plot_cells_trajectory.sh {input} {output} >& {log}'

rule plot_dynverse_cont_stratified_celltype:
    input:
        'output/hpbase/cont/seurat_annotated_landscaper.RData',
        'output/hpbase/cont_stratified/dynverse/{method}_celltype.RData'
    output:
        'plot/hpbase/cont_stratified/dynverse_{method}_celltype_clusters.png',
        'plot/hpbase/cont_stratified/dynverse_{method}_celltype_sample.png',
        'plot/hpbase/cont_stratified/dynverse_{method}_celltype_celltype.png'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_dynverse_cont_stratified_{method}.txt'
    log:
        'logs/plot_dynverse_cont_stratified_{method}.log'
    shell:
        'src/plot_dynverse.sh {input} {output} >& {log}'

rule plot_dynverse_DAPT_stratified_celltype:
    input:
        'output/hpbase/DAPT_stratified/seurat_annotated_landscaper.RData',
        'output/hpbase/DAPT_stratified/dynverse/{method}_celltype.RData'
    output:
        'plot/hpbase/DAPT_stratified/dynverse_{method}_celltype_clusters.png',
        'plot/hpbase/DAPT_stratified/dynverse_{method}_celltype_sample.png',
        'plot/hpbase/DAPT_stratified/dynverse_{method}_celltype_celltype.png'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_dynverse_DAPT_stratified_{method}.txt'
    log:
        'logs/plot_dynverse_DAPT_stratified_{method}.log'
    shell:
        'src/plot_dynverse.sh {input} {output} >& {log}'

#################################
# RNA Velocity
#################################
rule scv_pl_proportions_cont:
    input:
        'output/hpbase/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_proportions_cont_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_proportions_cont_dynamical.txt'
    log:
        'logs/scv_pl_proportions_cont_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_proportions_DAPT:
    input:
        'output/hpbase/DAPT/velocyto/DAPT_dynamical.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_proportions_DAPT_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_proportions_DAPT_dynamical.txt'
    log:
        'logs/scv_pl_proportions_DAPT_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream_cont:
    input:
        'output/hpbase/cont/velocyto/cont_{mode}.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_velocity_embedding_stream_cont_{mode}.png'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_cont_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_cont_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream_DAPT:
    input:
        'output/hpbase/DAPT/velocyto/DAPT_{mode}.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_velocity_embedding_stream_DAPT_{mode}.png'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_DAPT_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_DAPT_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_cont:
    input:
        'output/hpbase/cont/velocyto/cont_{mode}.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_velocity_embedding_cont_{mode}.png'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_cont_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_cont_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_DAPT:
    input:
        'output/hpbase/DAPT/velocyto/DAPT_{mode}.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_velocity_embedding_DAPT_{mode}.png'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_DAPT_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_DAPT_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_latenttime_cont:
    input:
        'output/hpbase/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_latenttime_cont_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_latenttime_cont_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_cont_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_latenttime_DAPT:
    input:
        'output/hpbase/DAPT/velocyto/DAPT_dynamical.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_latenttime_DAPT_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_latenttime_DAPT_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_DAPT_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_heatmap_cont:
    input:
        'output/hpbase/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_heatmap_cont_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_heatmap_cont_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_cont_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_heatmap_DAPT:
    input:
        'output/hpbase/DAPT/velocyto/DAPT_dynamical.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_heatmap_DAPT_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_heatmap_DAPT_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_DAPT_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers_cont:
    input:
        'output/hpbase/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_velocity_markers_cont_dynamical_1.png',
        'plot/hpbase/cont/scv_pl_velocity_markers_cont_dynamical_2.png',
        'plot/hpbase/cont/scv_pl_velocity_markers_cont_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_cont_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_cont_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers_DAPT:
    input:
        'output/hpbase/DAPT/velocyto/DAPT_dynamical.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_velocity_markers_DAPT_dynamical_1.png',
        'plot/hpbase/DAPT/scv_pl_velocity_markers_DAPT_dynamical_2.png',
        'plot/hpbase/DAPT/scv_pl_velocity_markers_DAPT_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_DAPT_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_DAPT_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'

rule scv_pl_paga_cont:
    input:
        'output/hpbase/cont/paga/cont.h5ad'
    output:
        'plot/hpbase/cont/scv_pl_paga_cont.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_paga_cont.txt'
    log:
        'logs/scv_pl_paga_cont.log'
    shell:
        'src/scv_pl_paga.sh {input} {output} >& {log}'

rule scv_pl_paga_DAPT:
    input:
        'output/hpbase/DAPT/paga/DAPT.h5ad'
    output:
        'plot/hpbase/DAPT/scv_pl_paga_DAPT.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scv_pl_paga_DAPT.txt'
    log:
        'logs/scv_pl_paga_DAPT.log'
    shell:
        'src/scv_pl_paga.sh {input} {output} >& {log}'
