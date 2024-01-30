import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        'output/hpbase/cont/seurat_celltype.RData',
        'output/hpbase/dapt/seurat_celltype.RData',
        expand('output/hpbase/{sample}/kana/{sample}.rds',
            sample=SAMPLES),
        'output/hpbase/integrated/kana/integrated.rds',
        expand('output/hpbase/{sample}/kana/stratification/{sample}.rds',
            sample=SAMPLES),
        'output/hpbase/cont/kana/cont.rds',
        'output/hpbase/DAPT/kana/DAPT.rds',
        'output/hpbase/24h/kana/24h.rds',
        'output/hpbase/36h/kana/36h.rds',
        'output/hpbase/48h/kana/48h.rds',
        'output/hpbase/72h/kana/72h.rds',
        'output/hpbase/96h/kana/96h.rds',
        expand('plot/hpbase/{sample}/dimplot_celltype.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dimplot_celltype.png',
        expand('plot/hpbase/{sample}/dimplot_germlayer.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dimplot_germlayer.png',
        'plot/hpbase/integrated/dimplot_celltype_splitby.png',
        expand('plot/hpbase/{sample}/dotplot_celltype.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/heatmap_celltype_cont.png',
        'plot/hpbase/integrated/heatmap_celltype_dapt.png',
        'plot/hpbase/integrated/heatmap_celltype_ratio.png',
        'plot/hpbase/integrated/heatmap_celltype_t_pval.png',
        'plot/hpbase/integrated/heatmap_celltype_wilcox_pval.png'

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

rule celltype_label_stratification:
    input:
        'output/hpbase/integrated/seurat_celltype.RData'
    output:
        'output/hpbase/cont/seurat_celltype.RData',
        'output/hpbase/dapt/seurat_celltype.RData'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/celltype_label_stratification.txt'
    log:
        'logs/celltype_label_stratification.log'
    shell:
        'src/celltype_label_stratification.sh {input} {output} >& {log}'

rule kana_sample:
    input:
        'output/hpbase/{sample}/seurat_celltype.RData'
    output:
        'output/hpbase/{sample}/kana/{sample}.rds'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/kana_{sample}.txt'
    log:
        'logs/kana_{sample}.log'
    shell:
        'src/kana.sh {input} {output} >& {log}'

rule kana_integrated:
    input:
        'output/hpbase/integrated/seurat_celltype.RData'
    output:
        'output/hpbase/integrated/kana/integrated.rds'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/kana_integrated.txt'
    log:
        'logs/kana_integrated.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_sample_stratification:
    input:
        'output/hpbase/integrated/kana/integrated.rds'
    output:
        'output/hpbase/{sample}/kana/stratification/{sample}.rds'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/kana_sample_stratification_{sample}.txt'
    log:
        'logs/kana_sample_stratification_{sample}.log'
    shell:
        'src/kana_sample_stratification.sh {input} {output} {wildcards.sample} >& {log}'

rule kana_experiment:
    input:
        'output/hpbase/integrated/kana/integrated.rds'
    output:
        'output/hpbase/cont/kana/cont.rds',
        'output/hpbase/DAPT/kana/DAPT.rds'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/kana_experiment.txt'
    log:
        'logs/kana_experiment.log'
    shell:
        'src/kana_experiment.sh {input} {output} >& {log}'

rule kana_time:
    input:
        'output/hpbase/integrated/kana/integrated.rds'
    output:
        'output/hpbase/24h/kana/24h.rds',
        'output/hpbase/36h/kana/36h.rds',
        'output/hpbase/48h/kana/48h.rds',
        'output/hpbase/72h/kana/72h.rds',
        'output/hpbase/96h/kana/96h.rds'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/kana_time.txt'
    log:
        'logs/kana_time.log'
    shell:
        'src/kana_time.sh {input} {output} >& {log}'

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

rule heatmap_celltype_integrated:
    input:
        'output/hpbase/integrated/seurat_celltype.RData'
    output:
        'plot/hpbase/integrated/heatmap_celltype_cont.png',
        'plot/hpbase/integrated/heatmap_celltype_dapt.png',
        'plot/hpbase/integrated/heatmap_celltype_ratio.png',
        'plot/hpbase/integrated/heatmap_celltype_t_pval.png',
        'plot/hpbase/integrated/heatmap_celltype_wilcox_pval.png'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/heatmap_celltype_integrated_hpbase_integrated.txt'
    log:
        'logs/heatmap_celltype_integrated_hpbase_integrated.log'
    shell:
        'src/heatmap_celltype_integrated.sh {input} {output} >& {log}'
