import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

rule all:
    input:
        'output/hpbase/integrated/seurat_scegot.h5ad',
        'output/hpbase/integrated/seurat_scegot_metadata.csv'

# #################################
# # Seurat => AnnData
# #################################
# rule seurat2anndata_integrated_scegot:
#     input:
#         'output/hpbase/integrated/seurat_annotated.RData'
#     output:
#         'output/hpbase/integrated/seurat_scegot.h5ad',
#         'output/hpbase/integrated/seurat_scegot_metadata.csv'
#     container:
#         'docker://koki/velocytor:20221015'
#     resources:
#         mem_mb=1000000
#     benchmark:
#         'benchmarks/seurat2anndata_integrated_scegot.txt'
#     log:
#         'logs/seurat2anndata_integrated_scegot.log'
#     shell:
#         'src/seurat2anndata_integrated_scegot.sh {input} {output} >& {log}'

#################################
# scEGOT
#################################
rule scegot:
    input:
        'output/hpbase/integrated/seurat_scegot.h5ad'
    output:
        'output/hpbase/integrated/scegot.h5ad'
    container:
        'docker://koki/scegot:20250515'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scegot.txt'
    log:
        'logs/scegot.log'
    shell:
        'src/scegot.sh {input} {output} >& {log}'

