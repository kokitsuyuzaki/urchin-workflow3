import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

DBS = ['hpbase', 'echinobase']

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        'data/hpbase/HpulGenome_v1_geneid.gtf',
        'data/echinobase/sp5_0_GCF_geneid.gtf',
        'data/geneid_to_genename.csv',
        'data/annotation.RData',
        'data/marker.RData',
        'data/go_bp_hpbase.RData',
        'data/go_mf_hpbase.RData',
        'data/go_cc_hpbase.RData',
        'data/go_bp_echinobase.RData',
        'data/go_mf_echinobase.RData',
        'data/go_cc_echinobase.RData'

#################################
# GTF files for RNA Velocity
#################################
rule gtf_hpbase:
    input:
        'data/hpbase/HpulGenome_v1.gtf'
    output:
        'data/hpbase/HpulGenome_v1_geneid.gtf'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/gtf_hpbase.txt'
    log:
        'logs/gtf_hpbase.log'
    shell:
        'src/gtf_geneid.sh {input} {output} >& {log}'

rule gtf_echinobase:
    input:
        'data/echinobase/sp5_0_GCF.gtf'
    output:
        'data/echinobase/sp5_0_GCF_geneid.gtf'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/gtf_echinobase.txt'
    log:
        'logs/gtf_echinobase.log'
    shell:
        'src/gtf_geneid.sh {input} {output} >& {log}'

#################################
# Corresponding Table for ID conversion
#################################
rule preprocess_geneid_to_genename:
    input:
        'data/hpbase/HpulGenome_v1_annot.xlsx'
    output:
        'data/geneid_to_genename.csv'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/preprocess_geneid_to_genename.txt'
    log:
        'logs/preprocess_geneid_to_genename.log'
    shell:
        'src/preprocess_geneid_to_genename.sh {input} {output} >& {log}'

#################################
# Gene Functional Annotation
#################################
rule preprocess_annotation:
    input:
        'data/hpbase/HpulGenome_v1_annot.xlsx'
    output:
        'data/annotation.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/preprocess_annotation.txt'
    log:
        'logs/preprocess_annotation.log'
    shell:
        'src/preprocess_annotation.sh {input} {output} >& {log}'

#################################
# Marker Gene table
#################################
rule preprocess_marker:
    input:
        'data/Shimoda/Genes_for_clustering.xlsx'
    output:
        'data/marker.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/preprocess_marker.txt'
    log:
        'logs/preprocess_marker.log'
    shell:
        'src/preprocess_marker.sh {input} {output} >& {log}'

#################################
# Gene Ontology
#################################
rule preprocess_go:
    input:
        'data/GeneGoTerms.txt'
    output:
        'data/go_bp_hpbase.RData',
        'data/go_mf_hpbase.RData',
        'data/go_cc_hpbase.RData',
        'data/go_bp_echinobase.RData',
        'data/go_mf_echinobase.RData',
        'data/go_cc_echinobase.RData'
    resources:
        mem_gb=100
    benchmark:
        'benchmarks/preprocess_go.txt'
    log:
        'logs/preprocess_go.log'
    shell:
        'src/preprocess_go.sh {input} {output} >& {log}'
