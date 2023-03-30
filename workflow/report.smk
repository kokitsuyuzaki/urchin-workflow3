import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
DBS = ['hpbase', 'echinobase']
GOS = ['bp', 'mf', 'cc']

container: 'docker://koki/urchin_workflow_seurat:20230111'

rule all:
    input:
        expand('plot/{db}/{sample}/scTGIF/{go}/index.html',
            db=DBS, sample=SAMPLES, go=GOS),
        expand('plot/{db}/integrated/scTGIF/{go}/index.html',
            db=DBS, go=GOS),
        expand('plot/{db}/cont/scTGIF/{go}/index.html',
            db=DBS, go=GOS),
        expand('plot/{db}/dapt/scTGIF/{go}/index.html',
            db=DBS, go=GOS)

#################################
# Functional Annotation
#################################
rule sctgif:
    input:
        'output/{db}/{sample}/seurat.RData',
        'data/go_{go}_{db}.RData'
    output:
        'plot/{db}/{sample}/scTGIF/{go}/index.html'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_{db}_{sample}_{go}.txt'
    log:
        'logs/sctgif_{db}_{sample}_{go}.log'
    shell:
        'src/sctgif.sh {wildcards.db} {wildcards.go} {input} {output} >& {log}'

rule sctgif_integrated:
    input:
        'output/{db}/integrated/seurat.RData',
        'data/go_{go}_{db}.RData'
    output:
        'plot/{db}/integrated/scTGIF/{go}/index.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_{db}_integrated_{go}.txt'
    log:
        'logs/sctgif_{db}_integrated_{go}.log'
    shell:
        'src/sctgif_integrated.sh {wildcards.db} {wildcards.go} {input} {output} >& {log}'

rule sctgif_cont:
    input:
        'output/{db}/cont/seurat.RData',
        'data/go_{go}_{db}.RData'
    output:
        'plot/{db}/cont/scTGIF/{go}/index.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_{db}_cont_{go}.txt'
    log:
        'logs/sctgif_{db}_cont_{go}.log'
    shell:
        'src/sctgif_cont.sh {wildcards.db} {wildcards.go} {input} {output} >& {log}'

rule sctgif_dapt:
    input:
        'output/{db}/dapt/seurat.RData',
        'data/go_{go}_{db}.RData'
    output:
        'plot/{db}/dapt/scTGIF/{go}/index.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sctgif_{db}_dapt_{go}.txt'
    log:
        'logs/sctgif_{db}_dapt_{go}.log'
    shell:
        'src/sctgif_dapt.sh {wildcards.db} {wildcards.go} {input} {output} >& {log}'
