import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
SAMPLES_STR = ['cont-24h_stratified', 'cont-36h_stratified', 'cont-48h_stratified', 'cont-72h_stratified', 'cont-96h_stratified', 'DAPT-24h_stratified', 'DAPT-36h_stratified', 'DAPT-48h_stratified', 'DAPT-72h_stratified', 'DAPT-96h_stratified']
CONDITIONS = ['cont', 'DAPT']
CONDITIONS_STR = ['cont_stratified', 'DAPT_stratified']
TIMES = ['24h', '36h', '48h', '72h', '96h']
TIMES_STR = ['24h_stratified', '36h_stratified', '48h_stratified', '72h_stratified', '96h_stratified']
CLUSTERS = list(range(41))
NEURONS_CLUSTERS = list(range(18))

container: 'docker://koki/urchin_workflow_seurat:20230616'

rule all:
    input:
        # Label transfer
        expand('output/hpbase/{sample}/seurat_annotated_lt.RData',
            sample=SAMPLES),
        # Kana
        'output/hpbase/integrated/kana/integrated.rds',
        'output/hpbase/integrated/kana/neurons.rds',
        expand('output/hpbase/integrated/kana/cluster{cl}.rds',
            cl=CLUSTERS),
        expand('output/hpbase/integrated/kana/neurons_cluster{ncl}.rds',
            ncl=NEURONS_CLUSTERS),
        expand('output/hpbase/{sample}/kana/{sample}.rds',
            sample=SAMPLES),
        expand('output/hpbase/{condition}/kana/{condition}.rds',
            condition=CONDITIONS),
        expand('output/hpbase/{time}/kana/{time}.rds',
            time=TIMES),
        expand('output/hpbase/{sample_str}/kana/{sample_str}.rds',
            sample_str=SAMPLES_STR),
        expand('output/hpbase/{condition_str}/kana/{condition_str}.rds',
            condition_str=CONDITIONS_STR),
        expand('output/hpbase/{time_str}/kana/{time_str}.rds',
            time_str=TIMES_STR),
        # Dimplot
        expand('plot/hpbase/{sample}/dimplot_celltype.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dimplot_celltype.png',
        'plot/hpbase/integrated/dimplot_celltype_splitby.png',
        expand('plot/hpbase/{sample}/dimplot_germlayer.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/dimplot_germlayer.png',
        'plot/hpbase/integrated/dimplot_neurons.png',
        'plot/hpbase/integrated/dimplot_neurons_splitby.png',
        'plot/hpbase/integrated/dimplot_cluster11.png',
        'plot/hpbase/integrated/dimplot_cluster11_splitby.png',
        # Dotplot
        expand('plot/hpbase/{sample}/dotplot_celltype.png',
            sample=SAMPLES),
        'plot/hpbase/integrated/heatmap_celltype_cont.png',
        'plot/hpbase/integrated/heatmap_celltype_DAPT.png',
        'plot/hpbase/integrated/heatmap_celltype_ratio.png',
        'plot/hpbase/integrated/heatmap_celltype_t_pval.png',
        'plot/hpbase/integrated/heatmap_celltype_wilcox_pval.png'

#################################
# Celltype/Germlayer Label
#################################
rule label_integrated:
    input:
        'output/hpbase/integrated/seurat.RData',
        'data/Shimoda/cluster_celltype.xlsx'
    output:
        'output/hpbase/integrated/seurat_annotated.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/label_integrated_hpbase_integrated.txt'
    log:
        'logs/label_integrated_hpbase_integrated.log'
    shell:
        'src/label_integrated.sh {input} {output} >& {log}'

# Stratification
rule label_integrated_neurons:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/seurat_annotated_neurons.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/label_integrated_hpbase_integrated_neurons.txt'
    log:
        'logs/label_integrated_hpbase_integrated_neurons.log'
    shell:
        'src/label_integrated_neurons.sh {input} {output} >& {log}'

rule label_integrated_neurons_clusters:
    input:
        'output/hpbase/integrated/seurat_annotated_neurons.RData'
    output:
        'output/hpbase/integrated/seurat_annotated_neurons_cluster{ncl}.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/label_integrated_hpbase_integrated_neurons_cluster{ncl}.txt'
    log:
        'logs/label_integrated_hpbase_integrated_neurons_cluster{ncl}.log'
    shell:
        'src/label_integrated_clusters.sh {input} {output} {wildcards.ncl} >& {log}'

rule label_integrated_clusters:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/seurat_annotated_cluster{cl}.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/label_integrated_hpbase_integrated_cluster{cl}.txt'
    log:
        'logs/label_integrated_hpbase_integrated_cluster{cl}.log'
    shell:
        'src/label_integrated_clusters.sh {input} {output} {wildcards.cl} >& {log}'

rule label_stratification:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        expand('output/hpbase/{sample}/seurat_annotated.RData',
            sample=SAMPLES),
        expand('output/hpbase/{condition}/seurat_annotated.RData',
            condition=CONDITIONS),
        expand('output/hpbase/{time}/seurat_annotated.RData',
            time=TIMES),
        expand('output/hpbase/{sample_str}/seurat_annotated.RData',
            sample_str=SAMPLES_STR),
        expand('output/hpbase/{condition_str}/seurat_annotated.RData',
            condition_str=CONDITIONS_STR),
        expand('output/hpbase/{time_str}/seurat_annotated.RData',
            time_str=TIMES_STR)
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/label_stratification.txt'
    log:
        'logs/label_stratification.log'
    shell:
        'src/label_stratification.sh >& {log}'

# For label transfer by urchin-integration-workflow
rule label_transfer:
    input:
        'data/geneid_to_genename.csv',
        'output/hpbase/{sample}/seurat_annotated.RData'
    output:
        'output/hpbase/{sample}/seurat_annotated_lt.RData'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/label_transfer_{sample}.txt'
    log:
        'logs/label_transfer_{sample}.log'
    shell:
        'src/label_transfer.sh {wildcards.sample} {input} {output} >& {log}'

# Kana
rule kana_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/kana/integrated.rds'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_integrated.txt'
    log:
        'logs/kana_integrated.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_integrated_neurons:
    input:
        'output/hpbase/integrated/seurat_annotated_neurons.RData'
    output:
        'output/hpbase/integrated/kana/neurons.rds'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_integrated_neurons.txt'
    log:
        'logs/kana_integrated_neurons.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_integrated_clusters:
    input:
        'output/hpbase/integrated/seurat_annotated_cluster{cl}.RData'
    output:
        'output/hpbase/integrated/kana/cluster{cl}.rds'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_integrated_cluster{cl}.txt'
    log:
        'logs/kana_integrated_cluster{cl}.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_integrated_neurons_clusters:
    input:
        'output/hpbase/integrated/seurat_annotated_neurons_cluster{cl}.RData'
    output:
        'output/hpbase/integrated/kana/neurons_cluster{cl}.rds'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_integrated_cluster{cl}.txt'
    log:
        'logs/kana_integrated_cluster{cl}.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_sample:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData'
    output:
        'output/hpbase/{sample}/kana/{sample}.rds'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_{sample}.txt'
    log:
        'logs/kana_{sample}.log'
    shell:
        'src/kana.sh {input} {output} >& {log}'

rule kana_condition:
    input:
        'output/hpbase/{condition}/seurat_annotated.RData'
    output:
        'output/hpbase/{condition}/kana/{condition}.rds'
    wildcard_constraints:
        condition='|'.join([re.escape(x) for x in CONDITIONS])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_condition_{condition}.txt'
    log:
        'logs/kana_condition_{condition}.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_time:
    input:
        'output/hpbase/{time}/seurat_annotated.RData'
    output:
        'output/hpbase/{time}/kana/{time}.rds'
    wildcard_constraints:
        time='|'.join([re.escape(x) for x in TIMES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_time_{time}.txt'
    log:
        'logs/kana_time_{time}.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_sample_stratified:
    input:
        'output/hpbase/{sample_str}/seurat_annotated.RData'
    output:
        'output/hpbase/{sample_str}/kana/{sample_str}.rds'
    wildcard_constraints:
        sample_str='|'.join([re.escape(x) for x in SAMPLES_STR])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_{sample_str}.txt'
    log:
        'logs/kana_{sample_str}.log'
    shell:
        'src/kana.sh {input} {output} >& {log}'

rule kana_condition_stratified:
    input:
        'output/hpbase/{condition_str}/seurat_annotated.RData'
    output:
        'output/hpbase/{condition_str}/kana/{condition_str}.rds'
    wildcard_constraints:
        condition_str='|'.join([re.escape(x) for x in CONDITIONS_STR])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_condition_stratified_{condition_str}.txt'
    log:
        'logs/kana_condition_stratified_{condition_str}.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

rule kana_time_stratified:
    input:
        'output/hpbase/{time_str}/seurat_annotated.RData'
    output:
        'output/hpbase/{time_str}/kana/{time_str}.rds'
    wildcard_constraints:
        time_str='|'.join([re.escape(x) for x in TIMES_STR])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/kana_time_stratified_{time_str}.txt'
    log:
        'logs/kana_time_stratified_{time_str}.log'
    shell:
        'src/kana_integrated.sh {input} {output} >& {log}'

#################################
# Dimplot
#################################
rule dimplot_celltype:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData'
    output:
        'plot/hpbase/{sample}/dimplot_celltype.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_celltype_hpbase_{sample}.txt'
    log:
        'logs/dimplot_celltype_hpbase_{sample}.log'
    shell:
        'src/dimplot_celltype.sh {input} {output} >& {log}'

rule dimplot_celltype_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'plot/hpbase/integrated/dimplot_celltype.png',
        'plot/hpbase/integrated/dimplot_celltype_splitby.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_celltype_integrated_hpbase_integrated.txt'
    log:
        'logs/dimplot_celltype_integrated_hpbase_integrated.log'
    shell:
        'src/dimplot_celltype_integrated.sh {input} {output} >& {log}'

rule dimplot_germlayer:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData'
    output:
        'plot/hpbase/{sample}/dimplot_germlayer.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_germlayer_hpbase_{sample}.txt'
    log:
        'logs/dimplot_germlayer_hpbase_{sample}.log'
    shell:
        'src/dimplot_germlayer.sh {input} {output} >& {log}'

rule dimplot_germlayer_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'plot/hpbase/integrated/dimplot_germlayer.png',
        'plot/hpbase/integrated/dimplot_germlayer_splitby.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_germlayer_integrated_hpbase_integrated.txt'
    log:
        'logs/dimplot_germlayer_integrated_hpbase_integrated.log'
    shell:
        'src/dimplot_germlayer_integrated.sh {input} {output} >& {log}'

rule dimplot_neurons_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated_neurons.RData'
    output:
        'plot/hpbase/integrated/dimplot_neurons.png',
        'plot/hpbase/integrated/dimplot_neurons_splitby.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_neurons_integrated_hpbase_integrated.txt'
    log:
        'logs/dimplot_neurons_integrated_hpbase_integrated.log'
    shell:
        'src/dimplot_neurons_integrated.sh {input} {output} >& {log}'

rule dimplot_cluster11_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated_cluster11.RData'
    output:
        'plot/hpbase/integrated/dimplot_cluster11.png',
        'plot/hpbase/integrated/dimplot_cluster11_splitby.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_cluster11_integrated_hpbase_integrated.txt'
    log:
        'logs/dimplot_cluster11_integrated_hpbase_integrated.log'
    shell:
        'src/dimplot_cluster11_integrated.sh {input} {output} >& {log}'

#################################
# Dotplot
#################################
# Figure 2
rule dotplot_celltype:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData'
    output:
        'plot/hpbase/{sample}/dotplot_celltype.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dotplot_celltype_hpbase_{sample}.txt'
    log:
        'logs/dotplot_celltype_hpbase_{sample}.log'
    shell:
        'src/dotplot_celltype.sh {input} {output} >& {log}'

# Figure 2
rule dotplot_celltype_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'plot/hpbase/integrated/dotplot_celltype.png',
        'plot/hpbase/integrated/dotplot_celltype_splitby.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dotplot_celltype_integrated_hpbase_integrated.txt'
    log:
        'logs/dotplot_celltype_integrated_hpbase_integrated.log'
    shell:
        'src/dotplot_celltype_integrated.sh {input} {output} >& {log}'

#################################
# Heatmap
#################################
# Figure 4
rule heatmap_celltype_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'plot/hpbase/integrated/heatmap_celltype_cont.png',
        'plot/hpbase/integrated/heatmap_celltype_DAPT.png',
        'plot/hpbase/integrated/heatmap_celltype_ratio.png',
        'plot/hpbase/integrated/heatmap_celltype_t_pval.png',
        'plot/hpbase/integrated/heatmap_celltype_wilcox_pval.png'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/heatmap_celltype_integrated_hpbase_integrated.txt'
    log:
        'logs/heatmap_celltype_integrated_hpbase_integrated.log'
    shell:
        'src/heatmap_celltype_integrated.sh {input} {output} >& {log}'
