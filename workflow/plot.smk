import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
DBS = ['hpbase', 'echinobase']
MODES = ['deterministic', 'stochastic', 'dynamical']
SAMPLE_PLOTS = ['elbowplot.png',
    'dimplot_cluster.png', 'featureplot_ncount_rna.png',
    'featureplot_nfeature_rna.png', 'featureplot_percent_mt.png',
    'featureplot_percent_rb.png',
    'ridgeplot_cellcycle.png', 'dimplot_cellcycle.png',
    'marker/FINISH_marker', 'marker/FINISH_cluster_marker',
    'featureplot_doublet.png', 
    'plot_cells_trajectory.png']
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

CONT_PLOTS = ['elbowplot.png', 'dimplot_cluster.png', 'dimplot_cluster_splitby.png', 'plot_cells_trajectory.png']

DAPT_PLOTS = CONT_PLOTS

rule all:
    input:
        # Basic Plots
        expand('plot/{db}/{sample}/{sample_plot}',
            sample=SAMPLES, db=DBS, sample_plot=SAMPLE_PLOTS),
        expand('plot/{db}/integrated/{integrated_plot}',
            db=DBS, integrated_plot=INTEGRATED_PLOTS),
        expand('plot/{db}/cont/{cont_plot}',
            db=DBS, cont_plot=CONT_PLOTS),
        expand('plot/{db}/dapt/{dapt_plot}',
            db=DBS, dapt_plot=DAPT_PLOTS),
        # Velocity Plots (Sample)
        expand('plot/{db}/{sample}/scv_pl_proportions_{sample}_dynamical.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/{sample}/scv_pl_velocity_embedding_stream_{sample}_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/{sample}/scv_pl_velocity_embedding_{sample}_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/{sample}/scv_pl_latenttime_{sample}_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/{sample}/scv_pl_heatmap_{sample}_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/{sample}/scv_pl_velocity_markers_{sample}_dynamical_1.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/{sample}/scv_pl_velocity_markers_{sample}_dynamical_2.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/{sample}/scv_pl_velocity_markers_{sample}_dynamical_3.png',
            db=DBS, sample=SAMPLES),
        # Velocity Plots (Integrated)
        expand('plot/{db}/integrated/scv_pl_proportions_integrated_dynamical.png',
            db=DBS, mode=MODES),
        expand('plot/{db}/integrated/scv_pl_velocity_embedding_stream_integrated_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/integrated/scv_pl_velocity_embedding_integrated_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/integrated/scv_pl_latenttime_integrated_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/integrated/scv_pl_heatmap_integrated_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/integrated/scv_pl_velocity_markers_integrated_dynamical_1.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/integrated/scv_pl_velocity_markers_integrated_dynamical_2.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/integrated/scv_pl_velocity_markers_integrated_dynamical_3.png',
            db=DBS, sample=SAMPLES),
        # Control
        expand('plot/{db}/cont/scv_pl_proportions_cont_dynamical.png',
            db=DBS, mode=MODES),
        expand('plot/{db}/cont/scv_pl_velocity_embedding_stream_cont_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/cont/scv_pl_velocity_embedding_cont_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/cont/scv_pl_latenttime_cont_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/cont/scv_pl_heatmap_cont_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/cont/scv_pl_velocity_markers_cont_dynamical_1.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/cont/scv_pl_velocity_markers_cont_dynamical_2.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/cont/scv_pl_velocity_markers_cont_dynamical_3.png',
            db=DBS, sample=SAMPLES),
        # DAPT
        expand('plot/{db}/dapt/scv_pl_proportions_dapt_dynamical.png',
            db=DBS, mode=MODES),
        expand('plot/{db}/dapt/scv_pl_velocity_embedding_stream_dapt_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/dapt/scv_pl_velocity_embedding_dapt_{mode}.png',
            db=DBS, sample=SAMPLES, mode=MODES),
        expand('plot/{db}/dapt/scv_pl_latenttime_dapt_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/dapt/scv_pl_heatmap_dapt_dynamical.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/dapt/scv_pl_velocity_markers_dapt_dynamical_1.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/dapt/scv_pl_velocity_markers_dapt_dynamical_2.png',
            db=DBS, sample=SAMPLES),
        expand('plot/{db}/dapt/scv_pl_velocity_markers_dapt_dynamical_3.png',
            db=DBS, sample=SAMPLES)

#################################
# Elbow Plot
#################################
rule elbowplot:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'plot/{db}/{sample}/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/elbowplot_{db}_{sample}.txt'
    log:
        'logs/elbowplot_{db}_{sample}.log'
    shell:
        'src/elbowplot.sh {input} {output} >& {log}'

rule elbowplot_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'plot/{db}/integrated/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/elbowplot_integrated_{db}_integrated.txt'
    log:
        'logs/elbowplot_integrated_{db}_integrated.log'
    shell:
        'src/elbowplot_integrated.sh {input} {output} >& {log}'

rule elbowplot_cont:
    input:
        'output/{db}/cont/seurat.RData'
    output:
        'plot/{db}/cont/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/elbowplot_cont_{db}_cont.txt'
    log:
        'logs/elbowplot_cont_{db}_cont.log'
    shell:
        'src/elbowplot_cont.sh {input} {output} >& {log}'

rule elbowplot_dapt:
    input:
        'output/{db}/dapt/seurat.RData'
    output:
        'plot/{db}/dapt/elbowplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/elbowplot_dapt_{db}_dapt.txt'
    log:
        'logs/elbowplot_dapt_{db}_dapt.log'
    shell:
        'src/elbowplot_dapt.sh {input} {output} >& {log}'

#################################
# Barplot
#################################
rule barplot_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'plot/{db}/integrated/barplot.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/barplot_integrated_{db}_integrated.txt'
    log:
        'logs/barplot_integrated_{db}_integrated.log'
    shell:
        'src/barplot_integrated.sh {input} {output} >& {log}'

#################################
# Cluster Label
#################################
rule dimplot_cluster:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'plot/{db}/{sample}/dimplot_cluster.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_cluster_{db}_{sample}.txt'
    log:
        'logs/dimplot_cluster_{db}_{sample}.log'
    shell:
        'src/dimplot_cluster.sh {input} {output} >& {log}'

rule dimplot_cluster_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'plot/{db}/integrated/dimplot_cluster.png',
        'plot/{db}/integrated/dimplot_cluster_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_cluster_integrated_{db}_integrated.txt'
    log:
        'logs/dimplot_cluster_integrated_{db}_integrated.log'
    shell:
        'src/dimplot_cluster_integrated.sh {input} {output} >& {log}'

rule dimplot_cluster_cont:
    input:
        'output/{db}/cont/seurat.RData'
    output:
        'plot/{db}/cont/dimplot_cluster.png',
        'plot/{db}/cont/dimplot_cluster_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_cluster_cont_{db}_cont.txt'
    log:
        'logs/dimplot_cluster_cont_{db}_cont.log'
    shell:
        'src/dimplot_cluster_cont.sh {input} {output} >& {log}'

rule dimplot_cluster_dapt:
    input:
        'output/{db}/dapt/seurat.RData'
    output:
        'plot/{db}/dapt/dimplot_cluster.png',
        'plot/{db}/dapt/dimplot_cluster_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_cluster_dapt_{db}_dapt.txt'
    log:
        'logs/dimplot_cluster_dapt_{db}_dapt.log'
    shell:
        'src/dimplot_cluster_dapt.sh {input} {output} >& {log}'

#################################
# Number of RNA counts
#################################
rule featureplot_ncount_rna:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'plot/{db}/{sample}/featureplot_ncount_rna.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_ncount_rna_{db}_{sample}.txt'
    log:
        'logs/featureplot_ncount_rna_{db}_{sample}.log'
    shell:
        'src/featureplot_ncount_rna.sh {input} {output} >& {log}'

rule featureplot_ncount_rna_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'plot/{db}/integrated/featureplot_ncount_rna.png',
        'plot/{db}/integrated/featureplot_ncount_rna_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_ncount_rna_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_ncount_rna_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_ncount_rna_integrated.sh {input} {output} >& {log}'

#################################
# Number of detected RNAs
#################################
rule featureplot_nfeature_rna:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'plot/{db}/{sample}/featureplot_nfeature_rna.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_nfeature_rna_{db}_{sample}.txt'
    log:
        'logs/featureplot_nfeature_rna_{db}_{sample}.log'
    shell:
        'src/featureplot_nfeature_rna.sh {input} {output} >& {log}'

rule featureplot_nfeature_rna_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'plot/{db}/integrated/featureplot_nfeature_rna.png',
        'plot/{db}/integrated/featureplot_nfeature_rna_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_nfeature_rna_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_nfeature_rna_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_nfeature_rna_integrated.sh {input} {output} >& {log}'

#################################
# Percentage of mitochondria genes' expression
#################################
rule featureplot_percent_mt:
    input:
        'output/{db}/{sample}/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/{db}/{sample}/featureplot_percent_mt.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_percent_mt_{db}_{sample}.txt'
    log:
        'logs/featureplot_percent_mt_{db}_{sample}.log'
    shell:
        'src/featureplot_percent_mt.sh {wildcards.db} {input} {output} >& {log}'

rule featureplot_percent_mt_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/{db}/integrated/featureplot_percent_mt.png',
        'plot/{db}/integrated/featureplot_percent_mt_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_percent_mt_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_percent_mt_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_percent_mt_integrated.sh {wildcards.db} {input} {output} >& {log}'

#################################
# Percentage of ribosome genes' expression
#################################
rule featureplot_percent_rb:
    input:
        'output/{db}/{sample}/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/{db}/{sample}/featureplot_percent_rb.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_percent_rb_{db}_{sample}.txt'
    log:
        'logs/featureplot_percent_rb_{db}_{sample}.log'
    shell:
        'src/featureplot_percent_rb.sh {wildcards.db} {input} {output} >& {log}'

rule featureplot_percent_rb_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/{db}/integrated/featureplot_percent_rb.png',
        'plot/{db}/integrated/featureplot_percent_rb_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_percent_rb_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_percent_rb_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_percent_rb_integrated.sh {wildcards.db} {input} {output} >& {log}'

#################################
# Cell cycle score
#################################
rule dimplot_cellcycle:
    input:
        'output/{db}/{sample}/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/{db}/{sample}/ridgeplot_cellcycle.png',
        'plot/{db}/{sample}/dimplot_cellcycle.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_cellcycle_{db}_{sample}.txt'
    log:
        'logs/dimplot_cellcycle_{db}_{sample}.log'
    shell:
        'src/dimplot_cellcycle.sh {wildcards.db} {input} {output} >& {log}'

rule dimplot_cellcycle_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'data/annotation.RData'
    output:
        'plot/{db}/integrated/ridgeplot_cellcycle.png',
        'plot/{db}/integrated/dimplot_cellcycle.png',
        'plot/{db}/integrated/dimplot_cellcycle_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_cellcycle_integrated_{db}_integrated.txt'
    log:
        'logs/dimplot_cellcycle_integrated_{db}_integrated.log'
    shell:
        'src/dimplot_cellcycle_integrated.sh {wildcards.db} {input} {output} >& {log}'

#################################
# Marker Genes (Prepared)
#################################
rule featureplot_marker:
    input:
        'output/{db}/{sample}/seurat.RData',
        'data/marker.RData'
    output:
        'plot/{db}/{sample}/marker/FINISH_marker'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_marker_{db}_{sample}.txt'
    log:
        'logs/featureplot_marker_{db}_{sample}.log'
    shell:
        'src/featureplot_marker.sh {wildcards.db} {input} {output} >& {log}'

rule featureplot_marker_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'data/marker.RData'
    output:
        'plot/{db}/integrated/marker/FINISH_marker'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_marker_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_marker_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_marker_integrated.sh {wildcards.db} {input} {output} >& {log}'

#################################
# Marker Genes (Detected in each Cluster)
#################################
rule featureplot_cluster_marker:
    input:
        'output/{db}/{sample}/seurat.RData',
        'output/{db}/{sample}/markers.xlsx'
    output:
        'plot/{db}/{sample}/marker/FINISH_cluster_marker'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_cluster_marker_{db}_{sample}.txt'
    log:
        'logs/featureplot_cluster_marker_{db}_{sample}.log'
    shell:
        'src/featureplot_cluster_marker.sh {input} {output} >& {log}'

rule featureplot_cluster_marker_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'output/{db}/integrated/markers.xlsx'
    output:
        'plot/{db}/integrated/marker/FINISH_cluster_marker'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_cluster_marker_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_cluster_marker_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_cluster_marker_integrated.sh {input} {output} >& {log}'

#################################
# Doublet Density
#################################
rule featureplot_doublet:
    input:
        'output/{db}/{sample}/seurat.RData',
        'output/{db}/{sample}/scdblfinder.RData'
    output:
        'plot/{db}/{sample}/featureplot_doublet.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_doublet_{db}_{sample}.txt'
    log:
        'logs/featureplot_doublet_{db}_{sample}.log'
    shell:
        'src/featureplot_doublet.sh {input} {output} >& {log}'

rule featureplot_doublet_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'output/{db}/integrated/scdblfinder.RData'
    output:
        'plot/{db}/integrated/featureplot_doublet.png',
        'plot/{db}/integrated/featureplot_doublet_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_doublet_integrated_{db}_integrated.txt'
    log:
        'logs/featureplot_doublet_integrated_{db}_integrated.log'
    shell:
        'src/featureplot_doublet_integrated.sh {input} {output} >& {log}'

#################################
# Trajectory Inference
#################################
rule plot_cells_trajectory:
    input:
        'output/{db}/{sample}/monocle3.RData'
    output:
        'plot/{db}/{sample}/plot_cells_trajectory.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/plot_cells_trajectory_{db}_{sample}.txt'
    log:
        'logs/plot_cells_trajectory_{db}_{sample}.log'
    shell:
        'src/plot_cells_trajectory.sh {input} {output} >& {log}'

rule plot_cells_traectory_integrated:
    input:
        'output/{db}/integrated/monocle3.RData'
    output:
        'plot/{db}/integrated/plot_cells_trajectory.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/plot_cells_trajectory_integrated_{db}_integrated.txt'
    log:
        'logs/plot_cells_trajectory_integrated_{db}_integrated.log'
    shell:
        'src/plot_cells_trajectory.sh {input} {output} >& {log}'

rule plot_cells_traectory_cont:
    input:
        'output/{db}/cont/monocle3.RData'
    output:
        'plot/{db}/cont/plot_cells_trajectory.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/plot_cells_trajectory_cont_{db}_cont.txt'
    log:
        'logs/plot_cells_trajectory_cont_{db}_cont.log'
    shell:
        'src/plot_cells_trajectory.sh {input} {output} >& {log}'

rule plot_cells_traectory_dapt:
    input:
        'output/{db}/dapt/monocle3.RData'
    output:
        'plot/{db}/dapt/plot_cells_trajectory.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/plot_cells_trajectory_dapt_{db}_dapt.txt'
    log:
        'logs/plot_cells_trajectory_dapt_{db}_dapt.log'
    shell:
        'src/plot_cells_trajectory.sh {input} {output} >& {log}'

#################################
# RNA Velocity
#################################
rule scv_pl_proportions:
    input:
        'output/{db}/{sample}/velocyto/{sample}_dynamical.h5ad'
    output:
        'plot/{db}/{sample}/scv_pl_proportions_{sample}_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    wildcard_constraints:
        db='|'.join([re.escape(x) for x in DBS]),
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        mode='|'.join([re.escape(x) for x in MODES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_proportions_{db}_{sample}_dynamical.txt'
    log:
        'logs/scv_pl_proportions_{db}_{sample}_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_proportions_integrated:
    input:
        'output/{db}/integrated/velocyto/integrated_dynamical.h5ad'
    output:
        'plot/{db}/integrated/scv_pl_proportions_integrated_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_proportions_{db}_integrated_dynamical.txt'
    log:
        'logs/scv_pl_proportions_{db}_integrated_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_proportions_cont:
    input:
        'output/{db}/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/{db}/cont/scv_pl_proportions_cont_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_proportions_{db}_cont_dynamical.txt'
    log:
        'logs/scv_pl_proportions_{db}_cont_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_proportions_dapt:
    input:
        'output/{db}/dapt/velocyto/dapt_dynamical.h5ad'
    output:
        'plot/{db}/dapt/scv_pl_proportions_dapt_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_proportions_{db}_dapt_dynamical.txt'
    log:
        'logs/scv_pl_proportions_{db}_dapt_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream:
    input:
        'output/{db}/{sample}/velocyto/{sample}_{mode}.h5ad'
    output:
        'plot/{db}/{sample}/scv_pl_velocity_embedding_stream_{sample}_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    wildcard_constraints:
        db='|'.join([re.escape(x) for x in DBS]),
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        mode='|'.join([re.escape(x) for x in MODES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_{db}_{sample}_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_{db}_{sample}_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream_integrated:
    input:
        'output/{db}/integrated/velocyto/integrated_{mode}.h5ad'
    output:
        'plot/{db}/integrated/scv_pl_velocity_embedding_stream_integrated_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_{db}_integrated_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_{db}_integrated_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream_cont:
    input:
        'output/{db}/cont/velocyto/cont_{mode}.h5ad'
    output:
        'plot/{db}/cont/scv_pl_velocity_embedding_stream_cont_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_{db}_cont_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_{db}_cont_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream_dapt:
    input:
        'output/{db}/dapt/velocyto/dapt_{mode}.h5ad'
    output:
        'plot/{db}/dapt/scv_pl_velocity_embedding_stream_dapt_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_{db}_dapt_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_{db}_dapt_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding:
    input:
        'output/{db}/{sample}/velocyto/{sample}_{mode}.h5ad'
    output:
        'plot/{db}/{sample}/scv_pl_velocity_embedding_{sample}_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    wildcard_constraints:
        db='|'.join([re.escape(x) for x in DBS]),
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        mode='|'.join([re.escape(x) for x in MODES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_{db}_{sample}_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_{db}_{sample}_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_integrated:
    input:
        'output/{db}/integrated/velocyto/integrated_{mode}.h5ad'
    output:
        'plot/{db}/integrated/scv_pl_velocity_embedding_integrated_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_{db}_integrated_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_{db}_integrated_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_cont:
    input:
        'output/{db}/cont/velocyto/cont_{mode}.h5ad'
    output:
        'plot/{db}/cont/scv_pl_velocity_embedding_cont_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_{db}_cont_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_{db}_cont_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_dapt:
    input:
        'output/{db}/dapt/velocyto/dapt_{mode}.h5ad'
    output:
        'plot/{db}/dapt/scv_pl_velocity_embedding_dapt_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_{db}_dapt_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_{db}_dapt_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_latenttime:
    input:
        'output/{db}/{sample}/velocyto/{sample}_dynamical.h5ad'
    output:
        'plot/{db}/{sample}/scv_pl_latenttime_{sample}_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    wildcard_constraints:
        db='|'.join([re.escape(x) for x in DBS]),
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        mode='|'.join([re.escape(x) for x in MODES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_latenttime_{db}_{sample}_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_{db}_{sample}_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_latenttime_integrated:
    input:
        'output/{db}/integrated/velocyto/integrated_dynamical.h5ad'
    output:
        'plot/{db}/integrated/scv_pl_latenttime_integrated_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_latenttime_{db}_integrated_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_{db}_integrated_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_latenttime_cont:
    input:
        'output/{db}/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/{db}/cont/scv_pl_latenttime_cont_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_latenttime_{db}_cont_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_{db}_cont_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_latenttime_dapt:
    input:
        'output/{db}/dapt/velocyto/dapt_dynamical.h5ad'
    output:
        'plot/{db}/dapt/scv_pl_latenttime_dapt_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_latenttime_{db}_dapt_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_{db}_dapt_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_heatmap:
    input:
        'output/{db}/{sample}/velocyto/{sample}_dynamical.h5ad'
    output:
        'plot/{db}/{sample}/scv_pl_heatmap_{sample}_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    wildcard_constraints:
        db='|'.join([re.escape(x) for x in DBS]),
        sample='|'.join([re.escape(x) for x in SAMPLES]),
        mode='|'.join([re.escape(x) for x in MODES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_heatmap_{db}_{sample}_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_{db}_{sample}_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_heatmap_integrated:
    input:
        'output/{db}/integrated/velocyto/integrated_dynamical.h5ad'
    output:
        'plot/{db}/integrated/scv_pl_heatmap_integrated_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_heatmap_{db}_integrated_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_{db}_integrated_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_heatmap_cont:
    input:
        'output/{db}/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/{db}/cont/scv_pl_heatmap_cont_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_heatmap_{db}_cont_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_{db}_cont_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_heatmap_dapt:
    input:
        'output/{db}/dapt/velocyto/dapt_dynamical.h5ad'
    output:
        'plot/{db}/dapt/scv_pl_heatmap_dapt_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_heatmap_{db}_dapt_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_{db}_dapt_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers:
    input:
        'output/{db}/{sample}/velocyto/{sample}_dynamical.h5ad'
    output:
        'plot/{db}/{sample}/scv_pl_velocity_markers_{sample}_dynamical_1.png',
        'plot/{db}/{sample}/scv_pl_velocity_markers_{sample}_dynamical_2.png',
        'plot/{db}/{sample}/scv_pl_velocity_markers_{sample}_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    wildcard_constraints:
        db='|'.join([re.escape(x) for x in DBS]),
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_{db}_{sample}_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_{db}_{sample}_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers_integrated:
    input:
        'output/{db}/integrated/velocyto/integrated_dynamical.h5ad'
    output:
        'plot/{db}/integrated/scv_pl_velocity_markers_integrated_dynamical_1.png',
        'plot/{db}/integrated/scv_pl_velocity_markers_integrated_dynamical_2.png',
        'plot/{db}/integrated/scv_pl_velocity_markers_integrated_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_{db}_integrated_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_{db}_integrated_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers_cont:
    input:
        'output/{db}/cont/velocyto/cont_dynamical.h5ad'
    output:
        'plot/{db}/cont/scv_pl_velocity_markers_cont_dynamical_1.png',
        'plot/{db}/cont/scv_pl_velocity_markers_cont_dynamical_2.png',
        'plot/{db}/cont/scv_pl_velocity_markers_cont_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_{db}_cont_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_{db}_cont_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers_dapt:
    input:
        'output/{db}/dapt/velocyto/dapt_dynamical.h5ad'
    output:
        'plot/{db}/dapt/scv_pl_velocity_markers_dapt_dynamical_1.png',
        'plot/{db}/dapt/scv_pl_velocity_markers_dapt_dynamical_2.png',
        'plot/{db}/dapt/scv_pl_velocity_markers_dapt_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_{db}_dapt_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_{db}_dapt_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'
