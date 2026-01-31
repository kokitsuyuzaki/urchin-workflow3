import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

METHODS = ["mst", "comp1", "angle"]

rule all:
    input:
        expand('output/hpbase/cont_stratified/dynverse/{method}_celltype.RData',
            method=METHODS),
        expand('output/hpbase/DAPT_stratified/dynverse/{method}_celltype.RData',
            method=METHODS),
        expand('output/hpbase/cont_stratified/dynverse/{method}_germlayer.RData',
            method=METHODS),
        expand('output/hpbase/DAPT_stratified/dynverse/{method}_germlayer.RData',
            method=METHODS),
        'output/hpbase/cont_stratified/monocle3.RData',
        'output/hpbase/DAPT_stratified/monocle3.RData',
        'output/hpbase/cont/paga/cont.h5ad',
        'output/hpbase/DAPT/paga/DAPT.h5ad'

#################################
# Dynverse
#################################
rule dynverse_cache:
    output:
        'cache_dir/{method}_FINISH'
    container:
        'docker://koki/dynverse:20251122'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/dynverse_cache_{method}.txt'
    log:
        'logs/dynverse_cache_{method}.log'
    shell:
        'src/dynverse_cache.sh {output} {wildcards.method} >& {log}'

rule start_id_cont_celltype:
    input:
        'output/hpbase/cont/seurat_annotated_landscaper.RData'
    output:
        'output/hpbase/cont_stratified/dynverse/start_id_celltype.RData'
    container:
        'docker://koki/dynverse:20251122'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/start_id_cont_celltype.txt'
    log:
        'logs/start_id_cont_celltype.log'
    shell:
        'src/start_id_celltype.sh {input} {output} >& {log}'

rule dynverse_cont_celltype:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/cont/seurat_annotated_landscaper.RData',
        'output/hpbase/cont_stratified/dynverse/start_id_celltype.RData'
    output:
        'output/hpbase/cont_stratified/dynverse/{method}_celltype.RData'
    container:
        'docker://koki/dynverse:20251122'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/dynverse_cont_{method}_celltype.txt'
    log:
        'logs/dynverse_cont_{method}_celltype.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

rule start_id_DAPT_celltype:
    input:
        'output/hpbase/DAPT/seurat_annotated_landscaper.RData'
    output:
        'output/hpbase/DAPT_stratified/dynverse/start_id_celltype.RData'
    container:
        'docker://koki/dynverse:20251122'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/start_id_DAPT_celltype.txt'
    log:
        'logs/start_id_DAPT_celltype.log'
    shell:
        'src/start_id_celltype.sh {input} {output} >& {log}'

rule dynverse_DAPT_celltype:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/DAPT/seurat_annotated_landscaper.RData',
        'output/hpbase/DAPT_stratified/dynverse/start_id_celltype.RData'
    output:
        'output/hpbase/DAPT_stratified/dynverse/{method}_celltype.RData'
    container:
        'docker://koki/dynverse:20251122'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/dynverse_DAPT_{method}_celltype.txt'
    log:
        'logs/dynverse_DAPT_{method}_celltype.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

#################################
# Monocle3
#################################
rule monocle3_cont:
    input:
        'output/hpbase/cont/seurat_annotated_landscaper.RData'
    output:
        'output/hpbase/cont_stratified/monocle3.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/monocle3_cont.txt'
    log:
        'logs/monocle3_cont.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'

rule monocle3_DAPT:
    input:
        'output/hpbase/DAPT/seurat_annotated_landscaper.RData'
    output:
        'output/hpbase/DAPT_stratified/monocle3.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000000
    benchmark:
        'benchmarks/monocle3_DAPT.txt'
    log:
        'logs/monocle3_DAPT.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'

#################################
# PAGA
#################################
rule paga_cont:
    input:
        'output/hpbase/cont/seurat_scegot_ectoderm.h5ad'
    output:
        'output/hpbase/cont/paga/cont.h5ad'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=100000000
    benchmark:
        'benchmarks/paga_cont.txt'
    log:
        'logs/paga_cont.log'
    shell:
        'src/paga.sh {input} {output} >& {log}'

rule paga_DAPT:
    input:
        'output/hpbase/DAPT/seurat_scegot_ectoderm.h5ad'
    output:
        'output/hpbase/DAPT/paga/DAPT.h5ad'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=100000000
    benchmark:
        'benchmarks/paga_DAPT.txt'
    log:
        'logs/paga_DAPT.log'
    shell:
        'src/paga.sh {input} {output} >& {log}'
