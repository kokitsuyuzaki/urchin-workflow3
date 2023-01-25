import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

rule all:
    input:
        'output/hpbase/integrated/tensorlycv/plot/bestrank_besttrial_barplot_FINISH',
        'output/hpbase/integrated/tensorlycv/plot/bestrank_besttrial_pairplot_FINISH'

#################################
# Functional Annotation
#################################
rule preprocess_tensorlycv:
    input:
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/urchin_tensor.npy'
    container:
        'docker://koki/urchin_workflow_seurat:20230111'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/preprocess_tensorlycv.txt'
    log:
        'logs/preprocess_tensorlycv.log'
    shell:
        'src/preprocess_tensorlycv.sh {input} {output} >& {log}'

rule tensorlycv:
    input:
        'output/hpbase/integrated/urchin_tensor.npy'
    output:
        'output/hpbase/integrated/tensorlycv/plot/bestrank_besttrial_barplot_FINISH',
        'output/hpbase/integrated/tensorlycv/plot/bestrank_besttrial_pairplot_FINISH'
    container:
        'docker://ghcr.io/kokitsuyuzaki/tensorlycv:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/tensorlycv.txt'
    log:
        'logs/tensorlycv.log'
    shell:
        'src/tensorlycv.sh {input} {output} >& {log}'
