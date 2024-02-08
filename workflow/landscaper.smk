import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ["integrated", "cont", "DAPT"]
PLOTFILES = ['ratio_group.png', 'Allstates.png', 'Freq_Prob_Energy.png', 'h.png', 'J.png', 'Basin.png', 'StatusNetwork_Subgraph.png', 'StatusNetwork_Subgraph_legend.png', 'StatusNetwork_Energy.png', 'StatusNetwork_Energy_legend.png', 'StatusNetwork_Ratio.png', 'StatusNetwork_Ratio_legend.png', 'StatusNetwork_State.png', 'StatusNetwork_State_legend.png', 'Landscape.png', 'discon_graph_1.png', 'discon_graph_2.png']

rule all:
    input:
        expand('plot/hpbase/{sample}/sbmfcv/deg/FINISH',
            sample=SAMPLES),
        'plot/hpbase/integrated/h.png',
        'plot/hpbase/integrated/J.png',
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
        'output/hpbase/integrated/seurat.RData'
    output:
        'output/hpbase/integrated/seurat.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/preprocess_sbmfcv.txt'
    log:
        'logs/preprocess_sbmfcv.log'
    shell:
        'src/preprocess_sbmfcv.sh {input} {output} >& {log}'

rule sbmfcv:
    input:
        'output/hpbase/integrated/seurat.tsv'
    output:
        'output/hpbase/integrated/sbmfcv/BIN_DATA.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/sbmfcv:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/sbmfcv.txt'
    log:
        'logs/sbmfcv.log'
    shell:
        'src/sbmfcv.sh {input} {output} >& {log}'

rule stratify:
    input:
        'output/hpbase/integrated/seurat.RData',
        'output/hpbase/integrated/sbmfcv/BIN_DATA.tsv'
    output:
        'output/hpbase/cont/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/DAPT/sbmfcv/BIN_DATA.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/stratify.txt'
    log:
        'logs/stratify.log'
    shell:
        'src/stratify.sh {input} {output} >& {log}'

rule deg_bindata:
    input:
        'output/hpbase/{sample}/seurat.RData',
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv'
    output:
        'output/hpbase/{sample}/sbmfcv/deg_bindata.xlsx'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/deg_bindata_{sample}.txt'
    log:
        'logs/deg_bindata_{sample}.log'
    shell:
        'src/deg_bindata_{wildcards.sample}.sh {input} {output} >& {log}'

rule featureplot_deg_bindata:
    input:
        'output/hpbase/{sample}/seurat.RData',
        'output/hpbase/{sample}/sbmfcv/deg_bindata.xlsx'
    output:
        'plot/hpbase/{sample}/sbmfcv/deg/FINISH'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_deg_bindata_{sample}.txt'
    log:
        'logs/featureplot_deg_bindata_{sample}.log'
    shell:
        'src/featureplot_deg_bindata_{wildcards.sample}.sh {input} {output} >& {log}'

rule group:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData'
    output:
        'output/hpbase/{sample}/group.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/group_{sample}.txt'
    log:
        'logs/group_{sample}.log'
    shell:
        'src/group_{wildcards.sample}.sh {input} {output} >& {log}'

rule landscaper_integrated:
    input:
        'output/hpbase/integrated/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/integrated/group.tsv'
    output:
        'plot/hpbase/integrated/Landscaper/plot/ratio_group.png',
        'plot/hpbase/integrated/Landscaper/plot/Allstates.png',
        'plot/hpbase/integrated/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/integrated/Landscaper/plot/h.png',
        'plot/hpbase/integrated/Landscaper/plot/J.png',
        'plot/hpbase/integrated/Landscaper/plot/Basin.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_State.png',
        'plot/hpbase/integrated/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/hpbase/integrated/Landscaper/plot/Landscape.png',
        'plot/hpbase/integrated/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/integrated/Landscaper/plot/discon_graph_2.png',
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/landscaper_integrated.txt'
    log:
        'logs/landscaper_integrated.log'
    shell:
        'src/landscaper_integrated.sh {input} {output} >& {log}'

rule landscaper_cont:
    input:
        'output/hpbase/cont/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/cont/group.tsv',
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv'
    output:
        'plot/hpbase/cont/Landscaper/plot/ratio_group.png',
        'plot/hpbase/cont/Landscaper/plot/Allstates.png',
        'plot/hpbase/cont/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/cont/Landscaper/plot/h.png',
        'plot/hpbase/cont/Landscaper/plot/J.png',
        'plot/hpbase/cont/Landscaper/plot/Basin.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_State.png',
        'plot/hpbase/cont/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/hpbase/cont/Landscaper/plot/Landscape.png',
        'plot/hpbase/cont/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/cont/Landscaper/plot/discon_graph_2.png'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/landscaper_cont.txt'
    log:
        'logs/landscaper_cont.log'
    shell:
        'src/landscaper_cont.sh {input} {output} >& {log}'

rule landscaper_DAPT:
    input:
        'output/hpbase/DAPT/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/DAPT/group.tsv',
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv'
    output:
        'plot/hpbase/DAPT/Landscaper/plot/ratio_group.png',
        'plot/hpbase/DAPT/Landscaper/plot/Allstates.png',
        'plot/hpbase/DAPT/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/DAPT/Landscaper/plot/h.png',
        'plot/hpbase/DAPT/Landscaper/plot/J.png',
        'plot/hpbase/DAPT/Landscaper/plot/Basin.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_State.png',
        'plot/hpbase/DAPT/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/hpbase/DAPT/Landscaper/plot/Landscape.png',
        'plot/hpbase/DAPT/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/DAPT/Landscaper/plot/discon_graph_2.png'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:main'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/landscaper_DAPT.txt'
    log:
        'logs/landscaper_DAPT.log'
    shell:
        'src/landscaper_DAPT.sh {input} {output} >& {log}'

rule plot_h:
    input:
        'plot/hpbase/cont/Landscaper/plot/h.png',
        'plot/hpbase/DAPT/Landscaper/plot/h.png'
    output:
        'plot/hpbase/integrated/h.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/plot_h.txt'
    log:
        'logs/plot_h.log'
    shell:
        'src/plot_h.sh {output} >& {log}'

rule plot_J:
    input:
        'plot/hpbase/cont/Landscaper/plot/J.png',
        'plot/hpbase/DAPT/Landscaper/plot/J.png'
    output:
        'plot/hpbase/integrated/J.png'
    container:
        'docker://koki/urchin_workflow_seurat:20230616'
    resources:
        mem_gb=500
    benchmark:
        'benchmarks/plot_J.txt'
    log:
        'logs/plot_J.log'
    shell:
        'src/plot_J.sh {output} >& {log}'

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
