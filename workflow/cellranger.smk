import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
DBS = ['hpbase', 'echinobase']

container: 'docker://litd/docker-cellranger:v7.0.0'

rule all:
    input:
        expand('output/{db}/aggr/outs/web_summary.html', db=DBS)

#################################
# Cell Ranger Mkref
#################################
rule cellranger_mkref_hpbase:
    input:
        'data/hpbase/trim.HpulGenome_v1.gtf',
        'data/hpbase/HpulGenome_v1_scaffold.fa'
    output:
        'data/hpbase/HpulGenome_v1/star/SA'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/cellranger_mkref_hpbase.txt'
    log:
        'logs/cellranger_mkref_hpbase.log'
    shell:
        'src/cellranger_mkref_hpbase.sh {input} {output} >& {log}'

rule cellranger_mkref_echinobase:
    input:
        'data/echinobase/trim.sp5_0_GCF.gtf',
        'data/echinobase/sp5_0_GCF_genomic.fa'
    output:
        'data/echinobase/sp5_0/star/SA'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/cellranger_mkref_echinobase.txt'
    log:
        'logs/cellranger_mkref_echinobase.log'
    shell:
        'src/cellranger_mkref_echinobase.sh {input} {output} >& {log}'

#################################
# Cell Ranger Count
#################################
rule cellranger_count_hpbase:
    input:
        'data/hpbase/HpulGenome_v1/star/SA',
        'data/Azenta/00_Rawdata/{sample}'
    output:
        'output/hpbase/{sample}/outs/web_summary.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/cellranger_count_hpbase_{sample}.txt'
    log:
        'logs/cellranger_count_hpbase_{sample}.log'
    shell:
        'src/cellranger_count_hpbase.sh {wildcards.sample} >& {log}'

rule cellranger_count_echinobase:
    input:
        'data/echinobase/sp5_0/star/SA',
        'data/Azenta/00_Rawdata/{sample}'
    output:
        'output/echinobase/{sample}/outs/web_summary.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/cellranger_count_echinobase_{sample}.txt'
    log:
        'logs/cellranger_count_echinobase_{sample}.log'
    shell:
        'src/cellranger_count_echinobase.sh {wildcards.sample} >& {log}'

#################################
# Cell Ranger Aggregate
#################################
def aggregate_sample(db):
    out = []
    for j in range(len(SAMPLES)):
        out.append('output/' + db[0] + '/' + SAMPLES[j] + '/outs/web_summary.html')
    return(out)

rule aggregation_csv:
    input:
        aggregate_sample
    output:
        'output/{db}/integrated/aggregation.csv'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/aggregation_csv_{db}.txt'
    log:
        'logs/aggregation_csv_{db}.log'
    shell:
        'src/aggregation_csv.sh {output} >& {log}'

rule cellranger_aggr:
    input:
        'output/{db}/integrated/aggregation.csv'
    output:
        'output/{db}/aggr/outs/web_summary.html'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/cellranger_aggr_{db}.txt'
    log:
        'logs/cellranger_aggr_{db}.log'
    shell:
        'src/cellranger_aggr.sh {wildcards.db} {input} >& {log}'
    