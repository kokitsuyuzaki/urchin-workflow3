import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

container: 'docker://litd/docker-cellranger:v7.0.0'

HPBASE_FILES = ['HpulGenome_v1_scaffold.fa', 'HpulGenome_v1_nucl.fa',
    'HpulGenome_v1_prot.fa', 'HpulTranscriptome.fa',
    'HpulTranscriptome_nucl.fa', 'HpulTranscriptome_prot.fa',
    'HpulGenome_v1_annot.xlsx', 'HpulGenome_v1.gff3']
ECHINOBASE_FILES = ['sp5_0_GCF_genomic.fa',
    'sp5_0_GCF.gff3', 'sp5_0_GCF_transcripts.fa']

rule all:
    input:
        expand('data/hpbase/{file}',
            file=HPBASE_FILES),
        expand('data/echinobase/{file}',
            file=ECHINOBASE_FILES),
        'data/GeneGoTerms.txt'

#################################
# Data download
#################################
rule download_hpbase:
    output:
        'data/hpbase/{file}'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/download_hpbase_{file}.txt'
    log:
        'logs/download_hpbase_{file}.log'
    shell:
        'src/download_hpbase.sh {wildcards.file} >& {log}'

rule download_echinobase:
    output:
        'data/echinobase/{file}'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/download_echinobase_{file}.txt'
    log:
        'logs/download_echinobase_{file}.log'
    shell:
        'src/download_echinobase.sh {wildcards.file} >& {log}'

rule download_go:
    output:
        'data/GeneGoTerms.txt'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/download_go.txt'
    log:
        'logs/download_go.log'
    shell:
        'src/download_go.sh {output} >& {log}'
