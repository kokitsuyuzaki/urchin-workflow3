import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']

rule all:
    input:
        'output/hpbase/an2/kana/integrated_an2.rds'

#################################
# Corresponding Table for ID conversion
#################################
rule preprocess_geneid_to_genename_an2:
    input:
        'data/hpbase/HpulGenome_v1_annot.xlsx'
    output:
        'data/geneid_to_genename_an2.csv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/preprocess_geneid_to_genename_an2.txt'
    log:
        'logs/preprocess_geneid_to_genename_an2.log'
    shell:
        'src/preprocess_geneid_to_genename_an2.sh {input} {output} >& {log}'

#################################
# Cell Ranger Mkref
#################################
rule cellranger_mkref_hpbase_an2:
    input:
        'data/hpbase/trim.HpulGenome_v1_an2.gtf',
        'data/hpbase/HpulGenome_v1_scaffold.fa'
    output:
        'data/hpbase/HpulGenome_v1_an2/star/SA'
    container:
        'docker://litd/docker-cellranger:v7.0.0'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/cellranger_mkref_hpbase_an2.txt'
    log:
        'logs/cellranger_mkref_hpbase_an2.log'
    shell:
        'src/cellranger_mkref_hpbase_an2.sh {input} {output} >& {log}'

#################################
# Cell Ranger Count
#################################
rule cellranger_count_hpbase_an2:
    input:
        'data/hpbase/HpulGenome_v1_an2/star/SA',
        'data/Azenta/00_Rawdata/{sample}'
    output:
        'output/hpbase/an2/{sample}/outs/web_summary.html'
    container:
        'docker://litd/docker-cellranger:v7.0.0'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/cellranger_count_hpbase_an2_{sample}.txt'
    log:
        'logs/cellranger_count_hpbase_an2_{sample}.log'
    shell:
        'src/cellranger_count_hpbase_an2.sh {wildcards.sample} >& {log}'

#################################
# Seurat
#################################
rule seurat_an2:
    input:
        'data/geneid_to_genename_an2.csv',
        'output/hpbase/an2/{sample}/outs/web_summary.html'
    output:
        'output/hpbase/an2/{sample}/seurat.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/seurat_an2_{sample}.txt'
    log:
        'logs/seurat_an2_{sample}.log'
    shell:
        'src/seurat_an2.sh {wildcards.sample} {output} >& {log}'

rule merge_seurat_integration_an2:
    input:
        expand('output/hpbase/an2/{sample}/seurat.RData', sample=SAMPLES),
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/an2/integrated/seurat_merged.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/merge_seurat_integration_an2.txt'
    log:
        'logs/merge_seurat_integration_an2.log'
    shell:
        'src/merge_seurat_integration_an2.sh {output} >& {log}'

rule kana_an2:
    input:
        'output/hpbase/an2/integrated/seurat_merged.RData'
    output:
        'output/hpbase/an2/kana/integrated_an2.rds'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_an2.txt'
    log:
        'logs/kana_an2.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'