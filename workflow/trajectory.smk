import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

METHODS = ["paga", "paga_tree", "comp1", "mst", "angle"]

rule all:
    input:
        expand('output/hpbase/integrated/dynverse/{method}_celltype.RData',
            method=METHODS),
        expand('output/hpbase/cont_stratified/dynverse/{method}_celltype.RData',
            method=METHODS),
        expand('output/hpbase/DAPT_stratified/dynverse/{method}_celltype.RData',
            method=METHODS),
        expand('output/hpbase/integrated/dynverse/{method}_germlayer.RData',
            method=METHODS),
        expand('output/hpbase/cont_stratified/dynverse/{method}_germlayer.RData',
            method=METHODS),
        expand('output/hpbase/DAPT_stratified/dynverse/{method}_germlayer.RData',
            method=METHODS),
        'output/hpbase/integrated/monocle3.RData',
        'output/hpbase/cont_stratified/monocle3.RData',
        'output/hpbase/DAPT_stratified/monocle3.RData'

#################################
# Dynverse
#################################
rule dynverse_cache:
    output:
        'cache_dir/{method}_FINISH'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_cache_{method}.txt'
    log:
        'logs/dynverse_cache_{method}.log'
    shell:
        'src/dynverse_cache.sh {output} {wildcards.method} >& {log}'

rule start_id_integrated_celltype:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/dynverse/start_id_celltype.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/start_id_integrated_celltype.txt'
    log:
        'logs/start_id_integrated_celltype.log'
    shell:
        'src/start_id_celltype.sh {input} {output} >& {log}'

rule dynverse_integrated_celltype:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/integrated/dynverse/start_id_celltype.RData'
    output:
        'output/hpbase/integrated/dynverse/{method}_celltype.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_integrated_{method}_celltype.txt'
    log:
        'logs/dynverse_integrated_{method}_celltype.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

rule start_id_cont_celltype:
    input:
        'output/hpbase/cont_stratified/seurat_annotated.RData'
    output:
        'output/hpbase/cont_stratified/dynverse/start_id_celltype.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/start_id_cont_celltype.txt'
    log:
        'logs/start_id_cont_celltype.log'
    shell:
        'src/start_id_celltype.sh {input} {output} >& {log}'

rule dynverse_cont_celltype:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/cont_stratified/seurat_annotated.RData',
        'output/hpbase/cont_stratified/dynverse/start_id_celltype.RData'
    output:
        'output/hpbase/cont_stratified/dynverse/{method}_celltype.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_cont_{method}_celltype.txt'
    log:
        'logs/dynverse_cont_{method}_celltype.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

rule start_id_DAPT_celltype:
    input:
        'output/hpbase/DAPT_stratified/seurat_annotated.RData'
    output:
        'output/hpbase/DAPT_stratified/dynverse/start_id_celltype.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/start_id_DAPT_celltype.txt'
    log:
        'logs/start_id_DAPT_celltype.log'
    shell:
        'src/start_id_celltype.sh {input} {output} >& {log}'

rule dynverse_DAPT_celltype:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/DAPT_stratified/seurat_annotated.RData',
        'output/hpbase/DAPT_stratified/dynverse/start_id_celltype.RData'
    output:
        'output/hpbase/DAPT_stratified/dynverse/{method}_celltype.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_DAPT_{method}_celltype.txt'
    log:
        'logs/dynverse_DAPT_{method}_celltype.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

rule start_id_integrated_germlayer:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/dynverse/start_id_germlayer.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/start_id_integrated_germlayer.txt'
    log:
        'logs/start_id_integrated_germlayer.log'
    shell:
        'src/start_id_germlayer.sh {input} {output} >& {log}'

rule dynverse_integrated_germlayer:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/integrated/dynverse/start_id_germlayer.RData'
    output:
        'output/hpbase/integrated/dynverse/{method}_germlayer.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_integrated_{method}_germlayer.txt'
    log:
        'logs/dynverse_integrated_{method}_germlayer.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

rule start_id_cont_germlayer:
    input:
        'output/hpbase/cont_stratified/seurat_annotated.RData'
    output:
        'output/hpbase/cont_stratified/dynverse/start_id_germlayer.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/start_id_cont_germlayer.txt'
    log:
        'logs/start_id_cont_germlayer.log'
    shell:
        'src/start_id_germlayer.sh {input} {output} >& {log}'

rule dynverse_cont_germlayer:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/cont_stratified/seurat_annotated.RData',
        'output/hpbase/cont_stratified/dynverse/start_id_germlayer.RData'
    output:
        'output/hpbase/cont_stratified/dynverse/{method}_germlayer.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_cont_{method}_germlayer.txt'
    log:
        'logs/dynverse_cont_{method}_germlayer.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

rule start_id_DAPT_germlayer:
    input:
        'output/hpbase/DAPT_stratified/seurat_annotated.RData'
    output:
        'output/hpbase/DAPT_stratified/dynverse/start_id_germlayer.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/start_id_DAPT_germlayer.txt'
    log:
        'logs/start_id_DAPT_germlayer.log'
    shell:
        'src/start_id_germlayer.sh {input} {output} >& {log}'

rule dynverse_DAPT_germlayer:
    input:
        'cache_dir/{method}_FINISH',
        'output/hpbase/DAPT_stratified/seurat_annotated.RData',
        'output/hpbase/DAPT_stratified/dynverse/start_id_germlayer.RData'
    output:
        'output/hpbase/DAPT_stratified/dynverse/{method}_germlayer.RData'
    container:
        'docker://koki/dynverse:20250203'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dynverse_DAPT_{method}_germlayer.txt'
    log:
        'logs/dynverse_DAPT_{method}_germlayer.log'
    shell:
        'src/dynverse.sh {input} {output} {wildcards.method} >& {log}'

#################################
# Monocle3
#################################
rule monocle3_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/monocle3.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/monocle3_integrated.txt'
    log:
        'logs/monocle3_integrated.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'

rule monocle3_cont:
    input:
        'output/hpbase/cont_stratified/seurat.RData'
    output:
        'output/hpbase/cont_stratified/monocle3.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/monocle3_cont.txt'
    log:
        'logs/monocle3_cont.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'

rule monocle3_DAPT:
    input:
        'output/hpbase/DAPT_stratified/seurat.RData'
    output:
        'output/hpbase/DAPT_stratified/monocle3.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/monocle3_DAPT.txt'
    log:
        'logs/monocle3_DAPT.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'
