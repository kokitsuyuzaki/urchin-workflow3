import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

container: 'docker://koki/urchin_workflow_bioconda:20220527'

rule all:
    input:
        'data/hpbase/trim.HpulGenome_v1.gtf',
        'data/echinobase/trim.sp5_0_GCF.gtf'

rule trim_gff_hpbase:
    input:
        'data/hpbase/HpulGenome_v1.gff3'
    output:
        'data/hpbase/trim.HpulGenome_v1.gff3'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/trim_gff_hpbase.txt'
    log:
        'logs/trim_gff_hpbase.log'
    shell:
        'src/trim_gff_hpbase.sh {input} {output} >& {log}'

rule gff_to_gtf:
    input:
        'data/hpbase/trim.HpulGenome_v1.gff3',
        'data/echinobase/sp5_0_GCF.gff3',
        'data/hpbase/HpulGenome_v1_scaffold.fa',
        'data/echinobase/sp5_0_GCF_genomic.fa'
    output:
        'data/hpbase/HpulGenome_v1.gtf',
        'data/echinobase/sp5_0_GCF.gtf'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/gff_to_gtf.txt'
    log:
        'logs/gff_to_gtf.log'
    shell:
        'src/gff_to_gtf.sh {input} {output} >& {log}'

rule trim_gtf:
    input:
        'data/hpbase/HpulGenome_v1.gtf',
        'data/echinobase/sp5_0_GCF.gtf'
    output:
        'data/hpbase/trim.HpulGenome_v1.gtf',
        'data/echinobase/trim.sp5_0_GCF.gtf'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/trim_gtf.txt'
    log:
        'logs/trim_gtf.log'
    shell:
        'src/trim_gtf.sh {input} {output} >& {log}'