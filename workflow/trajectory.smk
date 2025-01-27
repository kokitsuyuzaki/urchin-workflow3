import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
CONDITIONS = ['cont', 'DAPT']

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        'output/hpbase/integrated/monocle3.RData',
        expand('output/hpbase/{condition}/monocle3.RData',
            condition=CONDITIONS)

#################################
# Monocle3
#################################
rule monocle3_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/monocle3.RData'
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
        'output/hpbase/cont/seurat.RData'
    output:
        'output/hpbase/cont/monocle3.RData'
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
        'output/hpbase/DAPT/seurat.RData'
    output:
        'output/hpbase/DAPT/monocle3.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/monocle3_DAPT.txt'
    log:
        'logs/monocle3_DAPT.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'
