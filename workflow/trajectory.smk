import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
DBS = ['hpbase', 'echinobase']

container: 'docker://koki/urchin_workflow_seurat:20230111'

rule all:
    input:
        expand('output/{db}/{sample}/monocle3.RData',
            db=DBS, sample=SAMPLES),
        expand('output/{db}/integrated/monocle3.RData',
            db=DBS)

#################################
# Monocle3
#################################
rule monocle3:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'output/{db}/{sample}/monocle3.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/monocle3_{db}_{sample}.txt'
    log:
        'logs/monocle3_{db}_{sample}.log'
    shell:
        'src/monocle3.sh {input} {output} >& {log}'

rule monocle3_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'output/{db}/integrated/monocle3.RData'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/monocle3_{db}_integrated.txt'
    log:
        'logs/monocle3_{db}_integrated.log'
    shell:
        'src/monocle3_integrated.sh {input} {output} >& {log}'
