import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
DBS = ['hpbase', 'echinobase']

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        expand('output/{db}/{sample}/seurat_lt.RData',
            db=DBS, sample=SAMPLES),
        expand('output/{db}/cont/seurat.RData',
            db=DBS),
        expand('output/{db}/dapt/seurat.RData',
            db=DBS),
        expand('output/{db}/{sample}/markers.xlsx',
            db=DBS, sample=SAMPLES),
        expand('output/{db}/integrated/markers.xlsx',
            db=DBS),
        expand('output/{db}/cont/markers.xlsx',
            db=DBS),
        expand('output/{db}/dapt/markers.xlsx',
            db=DBS)

#################################
# Seurat
#################################
rule seurat:
    input:
        'data/geneid_to_genename.csv',
        'output/{db}/{sample}/outs/web_summary.html'
    output:
        'output/{db}/{sample}/seurat.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_{db}_{sample}.txt'
    log:
        'logs/seurat_{db}_{sample}.log'
    shell:
        'src/seurat.sh {wildcards.db} {wildcards.sample} {output} >& {log}'

rule seurat_for_labeltransfer:
    input:
        'data/geneid_to_genename.csv',
        'output/{db}/{sample}/outs/web_summary.html'
    output:
        'output/{db}/{sample}/seurat_lt.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_{db}_{sample}.txt'
    log:
        'logs/seurat_{db}_{sample}.log'
    shell:
        'src/seurat_lt.sh {wildcards.db} {wildcards.sample} {output} >& {log}'

def aggregate_sample(db):
    out = []
    for j in range(len(SAMPLES)):
        out.append('output/' + db[0] + '/' + SAMPLES[j] + '/seurat.RData')
    return(out)

rule seurat_integration:
    input:
        aggregate_sample
    output:
        'output/{db}/integrated/seurat.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_{db}_integration.txt'
    log:
        'logs/seurat_{db}_integration.log'
    shell:
        'src/seurat_integration.sh {wildcards.db} {output} >& {log}'

rule seurat_stratification:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'output/{db}/cont/seurat.RData',
        'output/{db}/dapt/seurat.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_{db}_integrated.txt'
    log:
        'logs/seurat_{db}_integrated.log'
    shell:
        'src/seurat_stratification.sh {input} {output} >& {log}'

rule seurat_findmarkers:
    input:
        'output/{db}/{sample}/seurat.RData'
    output:
        'output/{db}/{sample}/markers.xlsx'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_findmarkers_{db}_{sample}.txt'
    log:
        'logs/seurat_findmarkers_{db}_{sample}.log'
    shell:
        'src/seurat_findmarkers.sh {input} {output} >& {log}'

rule seurat_findconservedmarkers:
    input:
        'output/{db}/integrated/seurat.RData'
    output:
        'output/{db}/integrated/markers.xlsx'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_findconservedmarkers_{db}_integrated.txt'
    log:
        'logs/seurat_findconservedmarkers_{db}_integrated.log'
    shell:
        'src/seurat_findconservedmarkers.sh {input} {output} >& {log}'

rule seurat_findconservedmarkers_cont:
    input:
        'output/{db}/cont/seurat.RData'
    output:
        'output/{db}/cont/markers.xlsx'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_findconservedmarkers_{db}_cont.txt'
    log:
        'logs/seurat_findconservedmarkers_{db}_cont.log'
    shell:
        'src/seurat_findconservedmarkers_cont.sh {input} {output} >& {log}'

rule seurat_findconservedmarkers_dapt:
    input:
        'output/{db}/dapt/seurat.RData'
    output:
        'output/{db}/dapt/markers.xlsx'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_findconservedmarkers_{db}_dapt.txt'
    log:
        'logs/seurat_findconservedmarkers_{db}_dapt.log'
    shell:
        'src/seurat_findconservedmarkers_dapt.sh {input} {output} >& {log}'
