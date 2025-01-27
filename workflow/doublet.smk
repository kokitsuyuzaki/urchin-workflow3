import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        'output/hpbase/integrated/scdblfinder.RData'

#################################
# Elbow Plot
#################################

rule scdblfinder_integrated:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/scdblfinder.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scdblfinder_integrated.txt'
    log:
        'logs/scdblfinder_integrated.log'
    shell:
        'src/scdblfinder_integrated.sh {input} {output} >& {log}'
