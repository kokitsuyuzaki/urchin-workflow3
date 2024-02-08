import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

GOS = ['bp', 'mf', 'cc']

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        expand('plot/hpbase/integrated/scTGIF/{go}/index.html',
            go=GOS),
        expand('plot/hpbase/cont/scTGIF/{go}/index.html',
            go=GOS),
        expand('plot/hpbase/DAPT/scTGIF/{go}/index.html',
            go=GOS)

#################################
# Functional Annotation
#################################
rule sctgif_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/go_{go}_hpbase.RData'
    output:
        'plot/hpbase/integrated/scTGIF/{go}/index.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_integrated_{go}.txt'
    log:
        'logs/sctgif_integrated_{go}.log'
    shell:
        'src/sctgif_integrated.sh {wildcards.go} {input} {output} >& {log}'

rule sctgif_cont:
    input:
        'output/hpbase/cont/seurat.RData',
        'data/go_{go}_hpbase.RData'
    output:
        'plot/hpbase/cont/scTGIF/{go}/index.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_cont_{go}.txt'
    log:
        'logs/sctgif_cont_{go}.log'
    shell:
        'src/sctgif_integrated.sh {wildcards.go} {input} {output} >& {log}'

rule sctgif_DAPT:
    input:
        'output/hpbase/DAPT/seurat.RData',
        'data/go_{go}_hpbase.RData'
    output:
        'plot/hpbase/DAPT/scTGIF/{go}/index.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_DAPT_{go}.txt'
    log:
        'logs/sctgif_DAPT_{go}.log'
    shell:
        'src/sctgif_integrated.sh {wildcards.go} {input} {output} >& {log}'
