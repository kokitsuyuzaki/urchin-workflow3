import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

DBS = ['hpbase']
MODES = ['deterministic', 'stochastic', 'dynamical']

rule all:
    input:
        expand('plot/{db}/wo_2496h/scv_pl_proportions_wo_2496h_dynamical.png',
            db=DBS),
        expand('plot/{db}/wo_2496h/scv_pl_velocity_embedding_stream_wo_2496h_{mode}.png',
            db=DBS, mode=MODES),
        expand('plot/{db}/wo_2496h/scv_pl_velocity_embedding_wo_2496h_{mode}.png',
            db=DBS, mode=MODES),
        expand('plot/{db}/wo_2496h/scv_pl_latenttime_wo_2496h_dynamical.png',
            db=DBS),
        expand('plot/{db}/wo_2496h/scv_pl_heatmap_wo_2496h_dynamical.png',
            db=DBS),
        expand('plot/{db}/wo_2496h/scv_pl_velocity_markers_wo_2496h_dynamical_1.png',
            db=DBS),
        expand('plot/{db}/wo_2496h/scv_pl_velocity_markers_wo_2496h_dynamical_2.png',
            db=DBS),
        expand('plot/{db}/wo_2496h/scv_pl_velocity_markers_wo_2496h_dynamical_3.png',
            db=DBS)
        # expand('output/{db}/wo_2496h/velocyto/wo_2496h_{mode}.h5ad',
        #     db=DBS, mode=MODES)

rule stratify_wo_2496h:
    input:
        'output/{db}/integrated/seurat.RData',
    output:
        'output/{db}/wo_2496h/seurat.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/stratify_wo_2496h_{db}_wo_2496h.txt'
    log:
        'logs/stratify_wo_2496h_{db}_wo_2496h.log'
    shell:
        'src/stratify_wo_2496h.sh {input} {output} >& {log}'

#################################
# Aggregate Loom files
#################################

rule aggr_loom_wo_2496h:
    output:
        'output/{db}/aggr_wo_2496h/velocyto/aggr.loom'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/aggr_loom_wo_2496h_{db}.txt'
    log:
        'logs/aggr_loom_wo_2496h_{db}.log'
    shell:
        'src/aggr_loom_wo_2496h.sh {wildcards.db} {output} >& {log}'

#################################
# Seurat => AnnData
#################################

rule seurat2anndata_wo_2496h:
    input:
        'output/{db}/wo_2496h/seurat.RData'
    output:
        'output/{db}/wo_2496h/seurat.h5ad'
    container:
        'docker://koki/velocytor:20221015'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/seurat2anndata_{db}_wo_2496h.txt'
    log:
        'logs/seurat2anndata_{db}_wo_2496h.log'
    shell:
        'src/seurat2anndata_wo_2496h.sh {input} {output} >& {log}'

#################################
# Calculte RNA Velocity
#################################
rule scvelo_wo_2496h:
    input:
        'output/{db}/aggr_wo_2496h/velocyto/aggr.loom',
        'output/{db}/wo_2496h/seurat.h5ad'
    output:
        'output/{db}/wo_2496h/velocyto/wo_2496h_{mode}.h5ad'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scvelo_{db}_wo_2496h_{mode}.txt'
    log:
        'logs/scvelo_{db}_wo_2496h_{mode}.log'
    shell:
        'src/scvelo.sh {wildcards.mode} {input} {output} >& {log}'

#################################
# Plot
#################################

rule scv_pl_proportions_wo_2496h:
    input:
        'output/{db}/wo_2496h/velocyto/wo_2496h_dynamical.h5ad'
    output:
        'plot/{db}/wo_2496h/scv_pl_proportions_wo_2496h_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_proportions_{db}_wo_2496h_dynamical.txt'
    log:
        'logs/scv_pl_proportions_{db}_wo_2496h_dynamical.log'
    shell:
        'src/scv_pl_proportions.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_stream_wo_2496h:
    input:
        'output/{db}/wo_2496h/velocyto/wo_2496h_{mode}.h5ad'
    output:
        'plot/{db}/wo_2496h/scv_pl_velocity_embedding_stream_wo_2496h_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_stream_{db}_wo_2496h_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_stream_{db}_wo_2496h_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding_stream.sh {input} {output} >& {log}'

rule scv_pl_velocity_embedding_wo_2496h:
    input:
        'output/{db}/wo_2496h/velocyto/wo_2496h_{mode}.h5ad'
    output:
        'plot/{db}/wo_2496h/scv_pl_velocity_embedding_wo_2496h_{mode}.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_embedding_{db}_wo_2496h_{mode}.txt'
    log:
        'logs/scv_pl_velocity_embedding_{db}_wo_2496h_{mode}.log'
    shell:
        'src/scv_pl_velocity_embedding.sh {input} {output} >& {log}'

rule scv_pl_latenttime_wo_2496h:
    input:
        'output/{db}/wo_2496h/velocyto/wo_2496h_dynamical.h5ad'
    output:
        'plot/{db}/wo_2496h/scv_pl_latenttime_wo_2496h_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_latenttime_{db}_wo_2496h_dynamical.txt'
    log:
        'logs/scv_pl_latenttime_{db}_wo_2496h_dynamical.log'
    shell:
        'src/scv_pl_latenttime.sh {input} {output} >& {log}'

rule scv_pl_heatmap_wo_2496h:
    input:
        'output/{db}/wo_2496h/velocyto/wo_2496h_dynamical.h5ad'
    output:
        'plot/{db}/wo_2496h/scv_pl_heatmap_wo_2496h_dynamical.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_heatmap_{db}_wo_2496h_dynamical.txt'
    log:
        'logs/scv_pl_heatmap_{db}_wo_2496h_dynamical.log'
    shell:
        'src/scv_pl_heatmap.sh {input} {output} >& {log}'

rule scv_pl_velocity_markers_wo_2496h:
    input:
        'output/{db}/wo_2496h/velocyto/wo_2496h_dynamical.h5ad'
    output:
        'plot/{db}/wo_2496h/scv_pl_velocity_markers_wo_2496h_dynamical_1.png',
        'plot/{db}/wo_2496h/scv_pl_velocity_markers_wo_2496h_dynamical_2.png',
        'plot/{db}/wo_2496h/scv_pl_velocity_markers_wo_2496h_dynamical_3.png'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/scv_pl_velocity_markers_{db}_wo_2496h_dynamical.txt'
    log:
        'logs/scv_pl_velocity_markers_{db}_wo_2496h_dynamical.log'
    shell:
        'src/scv_pl_velocity_markers.sh {input} {output} >& {log}'