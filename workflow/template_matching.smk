import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

container: 'docker://koki/urchin_workflow_seurat:20230616'

TEMPLATE = ["Hp-Opn5L", "Hp-Tph", "Hp-Delta"]

rule all:
    input:
        expand('output/hpbase/integrated/template_matching_{template}.xlsx', template=TEMPLATE),
        expand('plot/hpbase/integrated/template_matching/{template}.png',
             template=TEMPLATE)

#################################
# Celltype Label
#################################
rule template_matching:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/template_matching_{template}.RData',
        'output/hpbase/integrated/template_matching_{template}.xlsx'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/template_matching_{template}.txt'
    log:
        'logs/template_matching_{template}.log'
    shell:
        'src/template_matching.sh {input} {output} {wildcards.template} >& {log}'

rule plot_template_matching:
    input:
        'output/hpbase/integrated/template_matching_{template}.RData'
    output:
        'plot/hpbase/integrated/template_matching/{template}.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/plot_template_matching_{template}.txt'
    log:
        'logs/plot_template_matching_{template}.log'
    shell:
        'src/plot_template_matching.sh {input} {output} {wildcards.template} >& {log}'
