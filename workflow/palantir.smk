import pandas as pd
from snakemake.utils import min_version

SAMPLES = ['cont', 'DAPT']
#################################
# Setting
#################################
min_version("6.5.3")

rule all:
    input:
        expand('plot/hpbase/{sample}/palantir.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/palantir_splitby.png',
            sample=SAMPLES)

#################################
# Palantir
#################################
rule palantir:
    input:
        'output/hpbase/{sample}/seurat_scegot_ectoderm.h5ad'
    output:
        'output/hpbase/{sample}/palantir_ectoderm.csv'
    container:
        'docker://koki/palantir:20251111'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/palantir_{sample}_ectoderm.txt'
    log:
        'logs/palantir_{sample}_ectoderm.log'
    shell:
        'src/palantir_{wildcards.sample}.sh {input} {output} >& {log}'

rule featureplot_palantir:
    input:
        'output/hpbase/{sample}/seurat_annotated_landscaper.RData',
        'output/hpbase/{sample}/palantir_ectoderm.csv'
    output:
        'plot/hpbase/{sample}/palantir.png',
        'plot/hpbase/{sample}/palantir_splitby.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_palantir_{sample}.txt'
    log:
        'logs/featureplot_palantir_{sample}.log'
    shell:
        'src/featureplot_palantir_{wildcards.sample}.sh {input} {output} >& {log}'
