import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
DBS = ['hpbase', 'echinobase']

container: 'docker://koki/urchin_workflow_seurat:20230111'

rule all:
    input:
        expand('output/{db}/{sample}/scdblfinder.RData',
            db=DBS, sample=SAMPLES),
        expand('output/{db}/integrated/scdblfinder.RData',
            db=DBS)

#################################
# Elbow Plot
#################################
rule scdblfinder:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'output/{db}/{sample}/scdblfinder.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/scdblfinder_{db}_{sample}.txt'
    log:
        'logs/scdblfinder_{db}_{sample}.log'
    shell:
        'src/scdblfinder.sh {input} {output} >& {log}'

rule scdblfinder_integrated:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'output/{db}/integrated/scdblfinder.RData'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/scdblfinder_{db}_integrated.txt'
    log:
        'logs/scdblfinder_{db}_integrated.log'
    shell:
        'src/scdblfinder_integrated.sh {input} {output} >& {log}'
