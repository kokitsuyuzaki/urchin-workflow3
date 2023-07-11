import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ["integrated", "cont", "dapt"]
PLOTFILES = ["Allstates.png", "Freq_Prob_Energy.png", "h.png", "J.png", "Basin.png", "StatusNetwork_Subgraph.png", "StatusNetwork_Energy.png", "Landscape.png", "discon_graph_1.png", "discon_graph_2.png"]

rule all:
    input:
        expand('plot/hpbase/{sample}/bindata/FINISH',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy_splitby.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/basin.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/basin_splitby.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/states.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/states_splitby.png',
            sample=SAMPLES)

#######################################
# Dimension Reduction & Binarization
#######################################
rule preprocess_sbmfcv:
    input:
        'output/hpbase/{sample}/seurat.RData'
    output:
        'output/hpbase/{sample}/seurat.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/preprocess_sbmfcv_{sample}.txt'
    log:
        'logs/preprocess_sbmfcv_{sample}.log'
    shell:
        'src/preprocess_sbmfcv_{wildcards.sample}.sh {input} {output} >& {log}'

rule sbmfcv:
    input:
        'output/hpbase/{sample}/seurat.tsv'
    output:
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/sbmfcv:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sbmfcv_{sample}.txt'
    log:
        'logs/sbmfcv_{sample}.log'
    shell:
        'src/sbmfcv.sh {input} {output} >& {log}'

rule landscaper:
    input:
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv'
    output:
        'plot/hpbase/{sample}/Landscaper/plot/Allstates.png',
        'plot/hpbase/{sample}/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/{sample}/Landscaper/plot/h.png',
        'plot/hpbase/{sample}/Landscaper/plot/J.png',
        'plot/hpbase/{sample}/Landscaper/plot/Basin.png',
        'plot/hpbase/{sample}/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/{sample}/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/{sample}/Landscaper/plot/Landscape.png',
        'plot/hpbase/{sample}/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/{sample}/Landscaper/plot/discon_graph_2.png'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/landscaper_{sample}.txt'
    log:
        'logs/landscaper_{sample}.log'
    shell:
        'src/landscaper.sh {input} {output} >& {log}'

rule dimplot_bindata:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/{sample}/bindata/FINISH'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/dimplot_bindata_{sample}.txt'
    log:
        'logs/dimplot_bindata_{sample}.log'
    shell:
        'src/dimplot_bindata_{wildcards.sample}.sh {output} >& {log}'

rule featureplot_energy:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/{sample}/energy.png',
        'plot/hpbase/{sample}/energy_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/featureplot_energy_{sample}.txt'
    log:
        'logs/featureplot_energy_{sample}.log'
    shell:
        'src/featureplot_energy_{wildcards.sample}.sh {output} >& {log}'

rule dimplot_basin:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/{sample}/basin.png',
        'plot/hpbase/{sample}/basin_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/dimplot_basin_{sample}.txt'
    log:
        'logs/dimplot_basin_{sample}.log'
    shell:
        'src/dimplot_basin_{wildcards.sample}.sh {output} >& {log}'

rule dimplot_states:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/{sample}/states.png',
        'plot/hpbase/{sample}/states_splitby.png',
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/dimplot_states_{sample}.txt'
    log:
        'logs/dimplot_states_{sample}.log'
    shell:
        'src/dimplot_states_{wildcards.sample}.sh {output} >& {log}'
