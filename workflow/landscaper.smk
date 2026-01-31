import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

rule all:
    input:
        'plot/hpbase/integrated/pca_landscaper.png',
        'plot/hpbase/integrated/featureplot_ncount_rna_landscaper.png'

#######################################
# Dimension Reduction & Binarization
#######################################
rule preprocess_landscaper:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/seurat_annotated_landscaper.RData',
        'output/hpbase/integrated/seurat.tsv',
        'output/hpbase/integrated/group.tsv',
        'output/hpbase/integrated_cov/group.tsv',
        'output/hpbase/cont/group.tsv',
        'output/hpbase/cont_cov/group.tsv',
        'output/hpbase/DAPT/group.tsv',
        'output/hpbase/DAPT_cov/group.tsv',
        'output/hpbase/integrated/cov.tsv',
        'output/hpbase/integrated_cov/cov.tsv',
        'output/hpbase/cont/cov.tsv',
        'output/hpbase/cont_cov/cov.tsv',
        'output/hpbase/DAPT/cov.tsv',
        'output/hpbase/DAPT_cov/cov.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_landscaper.txt'
    log:
        'logs/preprocess_landscaper.log'
    shell:
        'src/preprocess_landscaper.sh >& {log}'

rule pca_landscaper:
    input:
        'output/hpbase/integrated/seurat_annotated_landscaper.RData'
    output:
        'plot/hpbase/integrated/pca_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/pca_landscaper.txt'
    log:
        'logs/pca_landscaper.log'
    shell:
        'src/pca_landscaper.sh {input} {output} >& {log}'

rule featureplot_ncount_rna_landscaper:
    input:
        'output/hpbase/integrated/seurat_annotated_landscaper.RData'
    output:
        'plot/hpbase/integrated/featureplot_ncount_rna_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_ncount_rna_landscaper.txt'
    log:
        'logs/featureplot_ncount_rna_landscaper.log'
    shell:
        'src/featureplot_ncount_rna_landscaper.sh {input} {output} >& {log}'
