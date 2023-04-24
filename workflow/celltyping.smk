import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']

container: 'docker://koki/urchin_workflow_seurat:20230111'

rule all:
    input:
        expand('plot/hpbase/{sample}/dimplot_celltype.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dimplot_celltype.png',
        expand('plot/hpbase/{sample}/dimplot_germlayer.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dimplot_germlayer.png',
        'plot/hpbase/integrated/dimplot_celltype_splitby.png',
        expand('plot/hpbase/{sample}/dotplot_celltype.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dotplot_celltype.png',
        'plot/hpbase/integrated/dotplot_celltype_splitby.png'

#################################
# Celltype Label
#################################
rule celltype_label:
    input:
        'output/hpbase/{sample}/seurat.RData',
        'data/Shimoda/cluster_celltype.xlsx'
    output:
        'output/hpbase/{sample}/seurat_celltype.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/celltype_label_hpbase_{sample}.txt'
    log:
        'logs/celltype_label_hpbase_{sample}.log'
    shell:
        'src/celltype_label.sh {input} {output} >& {log}'

rule celltype_label_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/Shimoda/cluster_celltype.xlsx'
    output:
        'output/hpbase/integrated/seurat_celltype.RData'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/celltype_label_integrated_hpbase_integrated.txt'
    log:
        'logs/celltype_label_integrated_hpbase_integrated.log'
    shell:
        'src/celltype_label_integrated.sh {input} {output} >& {log}'

rule germlayer_label:
    input:
        'output/hpbase/{sample}/seurat.RData',
        'data/Shimoda/cluster_celltype.xlsx'
    output:
        'output/hpbase/{sample}/seurat_germlayer.RData'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/germlayer_label_hpbase_{sample}.txt'
    log:
        'logs/germlayer_label_hpbase_{sample}.log'
    shell:
        'src/germlayer_label.sh {input} {output} >& {log}'

rule germlayer_label_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/Shimoda/cluster_celltype.xlsx'
    output:
        'output/hpbase/integrated/seurat_germlayer.RData'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/germlayer_label_integrated_hpbase_integrated.txt'
    log:
        'logs/germlayer_label_integrated_hpbase_integrated.log'
    shell:
        'src/germlayer_label_integrated.sh {input} {output} >& {log}'

#################################
# Dimplot
#################################
rule dimplot_celltype:
    input:
        'output/hpbase/{sample}/seurat_celltype.RData'
    output:
        'plot/hpbase/{sample}/dimplot_celltype.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_celltype_hpbase_{sample}.txt'
    log:
        'logs/dimplot_celltype_hpbase_{sample}.log'
    shell:
        'src/dimplot_celltype.sh {input} {output} >& {log}'

rule dimplot_celltype_integrated:
    input:
        'output/hpbase/integrated/seurat_celltype.RData'
    output:
        'plot/hpbase/integrated/dimplot_celltype.png',
        'plot/hpbase/integrated/dimplot_celltype_splitby.png'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_celltype_integrated_hpbase_integrated.txt'
    log:
        'logs/dimplot_celltype_integrated_hpbase_integrated.log'
    shell:
        'src/dimplot_celltype_integrated.sh {input} {output} >& {log}'

rule dimplot_germlayer:
    input:
        'output/hpbase/{sample}/seurat_germlayer.RData'
    output:
        'plot/hpbase/{sample}/dimplot_germlayer.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_germlayer_hpbase_{sample}.txt'
    log:
        'logs/dimplot_germlayer_hpbase_{sample}.log'
    shell:
        'src/dimplot_celltype.sh {input} {output} >& {log}'

rule dimplot_germlayer_integrated:
    input:
        'output/hpbase/integrated/seurat_germlayer.RData'
    output:
        'plot/hpbase/integrated/dimplot_germlayer.png',
        'plot/hpbase/integrated/dimplot_germlayer_splitby.png'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dimplot_germlayer_integrated_hpbase_integrated.txt'
    log:
        'logs/dimplot_germlayer_integrated_hpbase_integrated.log'
    shell:
        'src/dimplot_celltype_integrated.sh {input} {output} >& {log}'

#################################
# Dotplot
#################################
rule dotplot_celltype:
    input:
        'output/hpbase/{sample}/seurat_celltype.RData'
    output:
        'plot/hpbase/{sample}/dotplot_celltype.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dotplot_celltype_hpbase_{sample}.txt'
    log:
        'logs/dotplot_celltype_hpbase_{sample}.log'
    shell:
        'src/dotplot_celltype.sh {input} {output} >& {log}'

rule dotplot_celltype_integrated:
    input:
        'output/hpbase/integrated/seurat_celltype.RData'
    output:
        'plot/hpbase/integrated/dotplot_celltype.png',
        'plot/hpbase/integrated/dotplot_celltype_splitby.png'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/dotplot_celltype_integrated_hpbase_integrated.txt'
    log:
        'logs/dotplot_celltype_integrated_hpbase_integrated.log'
    shell:
        'src/dotplot_celltype_integrated.sh {input} {output} >& {log}'