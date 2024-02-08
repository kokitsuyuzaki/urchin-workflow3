import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
CONDITIONS = ['cont', 'DAPT']
TIMES = ['24h', '36h', '48h', '72h', '96h']

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        expand('output/hpbase/{condition}/seurat.RData',
            condition=CONDITIONS),
        expand('output/hpbase/{time}/seurat.RData',
            time=TIMES),
        expand('output/echinobase/{sample}/seurat_lt.RData',
            sample=SAMPLES),
        'output/hpbase/integrated/markers.xlsx',
        'output/hpbase/cont_stratified/seurat.RData',
        'output/hpbase/DAPT_stratified/seurat.RData',
        'output/hpbase/24h_stratified/seurat.RData',
        'output/hpbase/36h_stratified/seurat.RData',
        'output/hpbase/48h_stratified/seurat.RData',
        'output/hpbase/72h_stratified/seurat.RData',
        'output/hpbase/96h_stratified/seurat.RData'

#################################
# Seurat
#################################
rule seurat:
    input:
        'data/geneid_to_genename.csv',
        'output/hpbase/{sample}/outs/web_summary.html'
    output:
        'output/hpbase/{sample}/seurat.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_{sample}.txt'
    log:
        'logs/seurat_{sample}.log'
    shell:
        'src/seurat.sh {wildcards.sample} {output} >& {log}'

######################################
# Seurat performed in each sample set
######################################
rule seurat_condition:
    input:
        expand('output/hpbase/{sample}/seurat.RData', sample=SAMPLES)
    output:
        'output/hpbase/{condition}/seurat.RData'
    resources:
        mem_gb=100
    wildcard_constraints:
        condition='|'.join([re.escape(x) for x in CONDITIONS])
    benchmark:
        'benchmarks/seurat_{condition}.txt'
    log:
        'logs/seurat_{condition}.log'
    shell:
        'src/seurat_condition.sh {wildcards.condition} {output} >& {log}'

rule seurat_time:
    input:
        expand('output/hpbase/{sample}/seurat.RData', sample=SAMPLES)
    output:
        'output/hpbase/{time}/seurat.RData'
    resources:
        mem_gb=100
    wildcard_constraints:
        time='|'.join([re.escape(x) for x in TIMES])
    benchmark:
        'benchmarks/seurat}_{time}.txt'
    log:
        'logs/seurat}_{time}.log'
    shell:
        'src/seurat_time.sh {wildcards.time} {output} >& {log}'

rule seurat_integration:
    input:
        expand('output/hpbase/{sample}/seurat.RData', sample=SAMPLES)
    output:
        'output/hpbase/integrated/seurat.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_integration.txt'
    log:
        'logs/seurat_integration.log'
    shell:
        'src/seurat_integration.sh {output} >& {log}'

rule seurat_for_labeltransfer:
    input:
        'data/geneid_to_genename.csv',
        'output/echinobase/{sample}/outs/web_summary.html'
    output:
        'output/echinobase/{sample}/seurat_lt.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_echinobase_{sample}.txt'
    log:
        'logs/seurat_echinobase_{sample}.log'
    shell:
        'src/seurat_lt.sh {wildcards.sample} {output} >& {log}'

rule seurat_findconservedmarkers_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/markers.xlsx'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_findconservedmarkers_integrated_integrated.txt'
    log:
        'logs/seurat_findconservedmarkers_integrated_integrated.log'
    shell:
        'src/seurat_findconservedmarkers.sh {input} {output} >& {log}'

################################################
# Stratification from integrated Seurat Object
################################################
rule seurat_stratification:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/cont_stratified/seurat.RData',
        'output/hpbase/DAPT_stratified/seurat.RData',
        'output/hpbase/24h_stratified/seurat.RData',
        'output/hpbase/36h_stratified/seurat.RData',
        'output/hpbase/48h_stratified/seurat.RData',
        'output/hpbase/72h_stratified/seurat.RData',
        'output/hpbase/96h_stratified/seurat.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/seurat_integrated.txt'
    log:
        'logs/seurat_integrated.log'
    shell:
        'src/seurat_stratification.sh {input} {output} >& {log}'
