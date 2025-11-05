import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ["integrated", "cont", "DAPT", "integrated_cov", "cont_cov", "DAPT_cov"]
SAMPLES2 = ["integrated", "cont_stratified", "DAPT_stratified"]
PLOTFILES = ['ratio_group.png', 'Allstates.png', 'Freq_Prob_Energy.png', 'h.png', 'J.png', 'Basin.png', 'StatusNetwork_Subgraph.png', 'StatusNetwork_Subgraph_legend.png', 'StatusNetwork_Energy.png', 'StatusNetwork_Energy_legend.png', 'StatusNetwork_Ratio.png', 'StatusNetwork_Ratio_legend.png', 'StatusNetwork_State.png', 'StatusNetwork_State_legend.png', 'Landscape.png', 'discon_graph_1.png', 'discon_graph_2.png']

rule all:
    input:
        'output/hpbase/integrated/sbmfcv/BIN_DATA.tsv',
        # expand('plot/hpbase/{sample}/sbmfcv/deg/FINISH',
        #     sample=SAMPLES),        
        expand('plot/hpbase/{sample}/bindata/FINISH',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/sce.RData',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy_hex.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy_rescaled.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy_rescaled_hex.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/energy_splitby.png',
            sample=SAMPLES),
        'plot/hpbase/cont_DAPT/energy_diff.png',
        'plot/hpbase/cont_DAPT_cov/energy_diff.png',
        expand('plot/hpbase/{sample2}/cellular_density.png',
            sample2=SAMPLES2),
        expand('plot/hpbase/{sample}/basin.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/basin_splitby.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/states.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/states_splitby.png',
            sample=SAMPLES),
        # 'plot/hpbase/cont/landscape_rescaled.png',
        # 'plot/hpbase/DAPT/landscape_rescaled.png',
        # 'plot/hpbase/cont_cov/landscape_rescaled.png',
        # 'plot/hpbase/DAPT_cov/landscape_rescaled.png',
        # Stable State Diagram
        'plot/hpbase/stable_states.png',
        # Transition Probability
        expand('plot/hpbase/{sample}/P_metropolis_rw.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_rw.png',
            sample=SAMPLES),
        # expand('plot/hpbase/{sample}/Landscaper/plot/discon_graph_1_new.png',
        #     sample=SAMPLES),
        # expand('plot/hpbase/{sample}/Landscaper/plot/discon_graph_2_new.png',
        #     sample=SAMPLES),
        # expand('plot/hpbase/{sample}/P_metropolis_emded.png',
        #     sample=SAMPLES),
        # expand('plot/hpbase/{sample}/P_glauber_emded.png',
        #     sample=SAMPLES),
        # Coarse-grained Vector Fields
        expand('plot/hpbase/{sample}/P_metropolis_cg.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_germlayer.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_cluster.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_sample.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_celltype.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_state.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_energy.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_germlayer.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_cluster.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_sample.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_celltype.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_state.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_energy.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_germlayer_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_cluster_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_sample_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_celltype_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_state_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_metropolis_cg_energy_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_germlayer_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_cluster_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_sample_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_celltype_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_state_streamline.png',
            sample=SAMPLES),
        expand('plot/hpbase/{sample}/P_glauber_cg_energy_streamline.png',
            sample=SAMPLES),

#######################################
# Dimension Reduction & Binarization
#######################################
rule preprocess_sbmfcv:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/seurat.tsv',
        'output/hpbase/integrated/group.tsv',
        'output/hpbase/integrated_cov/group.tsv',
        'output/hpbase/cont/group.tsv',
        'output/hpbase/cont_cov/group.tsv',
        'output/hpbase/DAPT/group.tsv',
        'output/hpbase/DAPT_cov/group.tsv',
        'output/hpbase/integrated/cov.tsv',
        'output/hpbase/integrated_cov/cov.tsv',
        'output/hpbase/cont/cov.tsv',
        'output/hpbase/cont_cov/cov.tsv',
        'output/hpbase/DAPT/cov.tsv',
        'output/hpbase/DAPT_cov/cov.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
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
        mem_mb=10000000
    benchmark:
        'benchmarks/sbmfcv.txt'
    log:
        'logs/sbmfcv.log'
    shell:
        'src/sbmfcv.sh {input} {output} >& {log}'

rule stratify:
    input:
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/integrated/sbmfcv/BIN_DATA.tsv'
    output:
        'output/hpbase/integrated_cov/seurat_annotated.RData',
        'output/hpbase/cont/seurat_annotated.RData',
        'output/hpbase/cont_cov/seurat_annotated.RData',
        'output/hpbase/DAPT/seurat_annotated.RData',
        'output/hpbase/DAPT_cov/seurat_annotated.RData',
        'output/hpbase/integrated_cov/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/cont/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/cont_cov/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/DAPT/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/DAPT_cov/sbmfcv/BIN_DATA.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/stratify.txt'
    log:
        'logs/stratify.log'
    shell:
        'src/stratify.sh {input} {output} >& {log}'

rule dimplot_bindata:
    input:
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv'
    output:
        'plot/hpbase/{sample}/bindata/FINISH'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_bindata_{sample}.txt'
    log:
        'logs/dimplot_bindata_{sample}.log'
    shell:
        'src/dimplot_bindata_{wildcards.sample}.sh {output} >& {log}'

rule deg_bindata:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData',
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv'
    output:
        'output/hpbase/{sample}/sbmfcv/deg_bindata.xlsx'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/deg_bindata_{sample}.txt'
    log:
        'logs/deg_bindata_{sample}.log'
    shell:
        'src/deg_bindata.sh {input} {output} >& {log}'

rule featureplot_deg_bindata:
    input:
        'output/hpbase/{sample}/seurat_annotated.RData',
        'output/hpbase/{sample}/sbmfcv/deg_bindata.xlsx'
    output:
        'plot/hpbase/{sample}/sbmfcv/deg/FINISH'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_gb=1000
    benchmark:
        'benchmarks/featureplot_deg_bindata_{sample}.txt'
    log:
        'logs/featureplot_deg_bindata_{sample}.log'
    shell:
        'src/featureplot_deg_bindata_{wildcards.sample}.sh {input} {output} >& {log}'

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
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv',
        'plot/hpbase/integrated/Landscaper/Allstates.tsv',
        'plot/hpbase/integrated/Landscaper/E.tsv',
        'plot/hpbase/integrated/Landscaper/Basin.tsv',
        'plot/hpbase/integrated/Landscaper/igraph.RData',
        'plot/hpbase/integrated/Landscaper/EnergyBarrier.tsv',
        'plot/hpbase/integrated/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
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
        'plot/hpbase/cont/Landscaper/plot/discon_graph_2.png',
        'plot/hpbase/cont/Landscaper/Basin.tsv',
        'plot/hpbase/cont/Landscaper/major_group.tsv',
        'plot/hpbase/cont/Landscaper/Allstates.tsv',
        'plot/hpbase/cont/Landscaper/E.tsv',
        'plot/hpbase/cont/Landscaper/igraph.RData',
        'plot/hpbase/cont/Landscaper/EnergyBarrier.tsv',
        'plot/hpbase/cont/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
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
        'plot/hpbase/DAPT/Landscaper/plot/discon_graph_2.png',
        'plot/hpbase/DAPT/Landscaper/Basin.tsv',
        'plot/hpbase/DAPT/Landscaper/major_group.tsv',
        'plot/hpbase/DAPT/Landscaper/Allstates.tsv',
        'plot/hpbase/DAPT/Landscaper/E.tsv',
        'plot/hpbase/DAPT/Landscaper/igraph.RData',
        'plot/hpbase/DAPT/Landscaper/EnergyBarrier.tsv',
        'plot/hpbase/DAPT/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_DAPT.txt'
    log:
        'logs/landscaper_DAPT.log'
    shell:
        'src/landscaper_DAPT.sh {input} {output} >& {log}'

rule landscaper_integrated_cov:
    input:
        'output/hpbase/integrated_cov/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/integrated_cov/group.tsv',
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv',
        'output/hpbase/integrated_cov/cov.tsv'
    output:
        'plot/hpbase/integrated_cov/Landscaper/plot/ratio_group.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/Allstates.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/h.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/J.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/Basin.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_State.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/Landscape.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/integrated_cov/Landscaper/plot/discon_graph_2.png',
        'plot/hpbase/integrated_cov/Landscaper/Allstates.tsv',
        'plot/hpbase/integrated_cov/Landscaper/E.tsv',
        'plot/hpbase/integrated_cov/Landscaper/Basin.tsv',
        'plot/hpbase/integrated_cov/Landscaper/igraph.RData',
        'plot/hpbase/integrated_cov/Landscaper/EnergyBarrier.tsv',
        'plot/hpbase/integrated_cov/Landscaper/Allstates_major_group.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_integrated_cov.txt'
    log:
        'logs/landscaper_integrated_cov.log'
    shell:
        'src/landscaper_integrated_cov.sh {input} {output} >& {log}'

rule landscaper_cont_cov:
    input:
        'output/hpbase/cont_cov/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/cont_cov/group.tsv',
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv',
        'output/hpbase/cont_cov/cov.tsv'
    output:
        'plot/hpbase/cont_cov/Landscaper/plot/ratio_group.png',
        'plot/hpbase/cont_cov/Landscaper/plot/Allstates.png',
        'plot/hpbase/cont_cov/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/cont_cov/Landscaper/plot/h.png',
        'plot/hpbase/cont_cov/Landscaper/plot/J.png',
        'plot/hpbase/cont_cov/Landscaper/plot/Basin.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_State.png',
        'plot/hpbase/cont_cov/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/hpbase/cont_cov/Landscaper/plot/Landscape.png',
        'plot/hpbase/cont_cov/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/cont_cov/Landscaper/plot/discon_graph_2.png',
        'plot/hpbase/cont_cov/Landscaper/h.tsv',
        'plot/hpbase/cont_cov/Landscaper/J.tsv',
        'plot/hpbase/cont_cov/Landscaper/g.txt',
        'plot/hpbase/cont_cov/Landscaper/Basin.tsv',
        'plot/hpbase/cont_cov/Landscaper/Allstates_major_group.tsv',
        'plot/hpbase/cont_cov/Landscaper/Allstates.tsv',
        'plot/hpbase/cont_cov/Landscaper/E.tsv',
        'plot/hpbase/cont_cov/Landscaper/igraph.RData',
        'plot/hpbase/cont_cov/Landscaper/EnergyBarrier.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_cont_cov.txt'
    log:
        'logs/landscaper_cont_cov.log'
    shell:
        'src/landscaper_cont_cov.sh {input} {output} >& {log}'

rule landscaper_DAPT_cov:
    input:
        'output/hpbase/DAPT_cov/sbmfcv/BIN_DATA.tsv',
        'output/hpbase/DAPT_cov/group.tsv',
        'plot/hpbase/integrated/Landscaper/Coordinate.tsv',
        'output/hpbase/DAPT_cov/cov.tsv'
    output:
        'plot/hpbase/DAPT_cov/Landscaper/plot/ratio_group.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/Allstates.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/h.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/J.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/Basin.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_State.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/Landscape.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/discon_graph_1.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/discon_graph_2.png',
        'plot/hpbase/DAPT_cov/Landscaper/h.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/J.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/g.txt',
        'plot/hpbase/DAPT_cov/Landscaper/Basin.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/Allstates_major_group.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/Allstates.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/E.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/igraph.RData',
        'plot/hpbase/DAPT_cov/Landscaper/EnergyBarrier.tsv'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/landscaper_DAPT_cov.txt'
    log:
        'logs/landscaper_DAPT_cov.log'
    shell:
        'src/landscaper_DAPT_cov.sh {input} {output} >& {log}'

rule transition_probability:
    input:
        'plot/hpbase/{sample}/Landscaper/E.tsv',
        'plot/hpbase/{sample}/Landscaper/Allstates.tsv'
    output:
        'plot/hpbase/{sample}/P_metropolis.tsv',
        'plot/hpbase/{sample}/P_glauber.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/transition_probability_{sample}.txt'
    log:
        'logs/transition_probability_{sample}.log'
    shell:
        'src/transition_probability.sh {input} {output} >& {log}'

rule coarse_grained_vector_fields:
    input:
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/{sample}/seurat_annotated.RData',
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv',
        'plot/hpbase/{sample}/Landscaper/Allstates.tsv',
        'plot/hpbase/{sample}/P_metropolis.tsv',
        'plot/hpbase/{sample}/P_glauber.tsv'
    output:
        'plot/hpbase/{sample}/P_metropolis_cg.RData',
        'plot/hpbase/{sample}/P_glauber_cg.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/coarse_grained_vector_fields_{sample}.txt'
    log:
        'logs/coarse_grained_vector_fields_{sample}.log'
    shell:
        'src/coarse_grained_vector_fields.sh {input} {output} >& {log}'

rule plot_coarse_grained_vector_fields:
    input:
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/{sample}/seurat_annotated.RData',
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv',
        'plot/hpbase/{sample}/Landscaper/Allstates.tsv',
        'plot/hpbase/{sample}/Landscaper/E.tsv',
        'plot/hpbase/{sample}/P_metropolis_cg.RData',
        'plot/hpbase/{sample}/P_glauber_cg.RData'
    output:
        'plot/hpbase/{sample}/P_metropolis_cg.png',
        'plot/hpbase/{sample}/P_metropolis_cg_germlayer.png',
        'plot/hpbase/{sample}/P_metropolis_cg_cluster.png',
        'plot/hpbase/{sample}/P_metropolis_cg_sample.png',
        'plot/hpbase/{sample}/P_metropolis_cg_celltype.png',
        'plot/hpbase/{sample}/P_metropolis_cg_state.png',
        'plot/hpbase/{sample}/P_metropolis_cg_energy.png',
        'plot/hpbase/{sample}/P_glauber_cg.png',
        'plot/hpbase/{sample}/P_glauber_cg_germlayer.png',
        'plot/hpbase/{sample}/P_glauber_cg_cluster.png',
        'plot/hpbase/{sample}/P_glauber_cg_sample.png',
        'plot/hpbase/{sample}/P_glauber_cg_celltype.png',
        'plot/hpbase/{sample}/P_glauber_cg_state.png',
        'plot/hpbase/{sample}/P_glauber_cg_energy.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_coarse_grained_vector_fields_{sample}.txt'
    log:
        'logs/plot_coarse_grained_vector_fields_{sample}.log'
    shell:
        'src/plot_coarse_grained_vector_fields.sh {input} {output} >& {log}'

rule plot_coarse_grained_vector_fields_streamline:
    input:
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/{sample}/seurat_annotated.RData',
        'output/hpbase/{sample}/sbmfcv/BIN_DATA.tsv',
        'plot/hpbase/{sample}/Landscaper/Allstates.tsv',
        'plot/hpbase/{sample}/Landscaper/E.tsv',
        'plot/hpbase/{sample}/P_metropolis_cg.RData',
        'plot/hpbase/{sample}/P_glauber_cg.RData'
    output:
        'plot/hpbase/{sample}/P_metropolis_cg_streamline.png',
        'plot/hpbase/{sample}/P_metropolis_cg_germlayer_streamline.png',
        'plot/hpbase/{sample}/P_metropolis_cg_cluster_streamline.png',
        'plot/hpbase/{sample}/P_metropolis_cg_sample_streamline.png',
        'plot/hpbase/{sample}/P_metropolis_cg_celltype_streamline.png',
        'plot/hpbase/{sample}/P_metropolis_cg_state_streamline.png',
        'plot/hpbase/{sample}/P_metropolis_cg_energy_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_germlayer_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_cluster_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_sample_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_celltype_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_state_streamline.png',
        'plot/hpbase/{sample}/P_glauber_cg_energy_streamline.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_coarse_grained_vector_fields_streamline_{sample}.txt'
    log:
        'logs/plot_coarse_grained_vector_fields_streamline_{sample}.log'
    shell:
        'src/plot_coarse_grained_vector_fields_streamline.sh {input} {output} >& {log}'

rule random_walk:
    input:
        'plot/hpbase/{sample}/P_metropolis.tsv',
        'plot/hpbase/{sample}/P_glauber.tsv'
    output:
        'plot/hpbase/{sample}/P_metropolis_rw.tsv',
        'plot/hpbase/{sample}/P_glauber_rw.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/random_walk_{sample}.txt'
    log:
        'logs/random_walk_{sample}.log'
    shell:
        'src/random_walk.sh {input} {output} >& {log}'

rule plot_random_walk_metropolis:
    input:
        'plot/hpbase/integrated/Landscaper/major_group.tsv',
        'plot/hpbase/{sample}/P_metropolis_rw.tsv'
    output:
        'plot/hpbase/{sample}/P_metropolis_rw.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_metropolis_{sample}.txt'
    log:
        'logs/plot_random_walk_metropolis_{sample}.log'
    shell:
        'src/plot_random_walk.sh {input} {output} >& {log}'

rule plot_random_walk_glauber:
    input:
        'plot/hpbase/integrated/Landscaper/major_group.tsv',
        'plot/hpbase/{sample}/P_glauber_rw.tsv'
    output:
        'plot/hpbase/{sample}/P_glauber_rw.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_glauber_{sample}.txt'
    log:
        'logs/plot_random_walk_glauber_{sample}.log'
    shell:
        'src/plot_random_walk.sh {input} {output} >& {log}'

rule graph_embedding:
    input:
        'plot/hpbase/{sample}/Landscaper/Allstates.tsv',
        'plot/hpbase/{sample}/P_metropolis.tsv',
        'plot/hpbase/{sample}/P_glauber.tsv'
    output:
        'plot/hpbase/{sample}/P_metropolis_emded.RData',
        'plot/hpbase/{sample}/P_glauber_emded.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/graph_embedding_{sample}.txt'
    log:
        'logs/graph_embedding_{sample}.log'
    shell:
        'src/graph_embedding.sh {input} {output} >& {log}'

rule plot_graph_embedding:
    input:
        'plot/hpbase/{sample}/Landscaper/E.tsv',
        'plot/hpbase/{sample}/Landscaper/Basin.tsv',
        'plot/hpbase/{sample}/Landscaper/igraph.RData',
        'plot/hpbase/{sample}/P_metropolis_emded.RData',
        'plot/hpbase/{sample}/P_glauber_emded.RData'
    output:
        'plot/hpbase/{sample}/P_metropolis_emded.png',
        'plot/hpbase/{sample}/P_glauber_emded.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_graph_embedding_{sample}.txt'
    log:
        'logs/plot_graph_embedding_{sample}.log'
    shell:
        'src/plot_graph_embedding.sh {input} {output} >& {log}'

rule tipping_state:
    input:
        'plot/hpbase/{sample}/Landscaper/Allstates.tsv',
        'plot/hpbase/{sample}/Landscaper/E.tsv',
        'plot/hpbase/{sample}/Landscaper/Basin.tsv'
    output:
        'plot/hpbase/{sample}/Landscaper/TippingState.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/tipping_state_{sample}.txt'
    log:
        'logs/tipping_state_{sample}.log'
    shell:
        'src/tipping_state.sh {input} {output} >& {log}'

rule dendrogram:
    input:
        'plot/hpbase/{sample}/Landscaper/E.tsv',
        'plot/hpbase/{sample}/Landscaper/Basin.tsv',
        'plot/hpbase/{sample}/Landscaper/EnergyBarrier.tsv',
        'plot/hpbase/{sample}/Landscaper/Allstates_major_group.tsv',
        'plot/hpbase/{sample}/Landscaper/TippingState.tsv'
    output:
        'plot/hpbase/{sample}/Landscaper/dendrogram_new.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dendrogram_{sample}.txt'
    log:
        'logs/dendrogram_{sample}.log'
    shell:
        'src/dendrogram.sh {input} {output} >& {log}'

rule plot_discon_graph:
    input:
        'plot/hpbase/{sample}/Landscaper/dendrogram_new.RData'
    output:
        'plot/hpbase/{sample}/Landscaper/plot/discon_graph_1_new.png',
        'plot/hpbase/{sample}/Landscaper/plot/discon_graph_2_new.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_discon_graph_{sample}.txt'
    log:
        'logs/plot_discon_graph_{sample}.log'
    shell:
        'src/plot_discon_graph.sh {input} {output} >& {log}'

rule plot_h:
    input:
        'plot/hpbase/cont/Landscaper/plot/h.png',
        'plot/hpbase/DAPT/Landscaper/plot/h.png',
        'plot/hpbase/cont_cov/Landscaper/plot/h.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/h.png'
    output:
        'plot/hpbase/integrated/h.png',
        'plot/hpbase/integrated_cov/h.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_h.txt'
    log:
        'logs/plot_h.log'
    shell:
        'src/plot_h.sh {output} >& {log}'

rule plot_J:
    input:
        'plot/hpbase/cont/Landscaper/plot/J.png',
        'plot/hpbase/DAPT/Landscaper/plot/J.png',
        'plot/hpbase/cont_cov/Landscaper/plot/J.png',
        'plot/hpbase/DAPT_cov/Landscaper/plot/J.png'
    output:
        'plot/hpbase/integrated/J.png',
        'plot/hpbase/integrated_cov/J.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_J.txt'
    log:
        'logs/plot_J.log'
    shell:
        'src/plot_J.sh {output} >& {log}'

rule featureplot_energy:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/{sample}/sce.RData',
        'plot/hpbase/{sample}/energy.png',
        'plot/hpbase/{sample}/energy_hex.png',
        'plot/hpbase/{sample}/energy_rescaled.png',
        'plot/hpbase/{sample}/energy_rescaled_hex.png',
        'plot/hpbase/{sample}/energy_splitby.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_{sample}.txt'
    log:
        'logs/featureplot_energy_{sample}.log'
    shell:
        'src/featureplot_energy_{wildcards.sample}.sh {output} >& {log}'

rule featureplot_energy_cont_dapt:
    input:
        'plot/hpbase/integrated/sce.RData',
        expand('plot/hpbase/cont/Landscaper/plot/{p}', p=PLOTFILES),
        expand('plot/hpbase/DAPT/Landscaper/plot/{p}', p=PLOTFILES)
    output:
        'plot/hpbase/cont_DAPT/energy_diff.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_cont_DAPT.txt'
    log:
        'logs/featureplot_energy_cont_DAPT.log'
    shell:
        'src/featureplot_energy_cont_DAPT.sh >& {log}'

rule featureplot_energy_cont_dapt_cov:
    input:
        'plot/hpbase/integrated_cov/sce.RData',
        expand('plot/hpbase/cont_cov/Landscaper/plot/{p}', p=PLOTFILES),
        expand('plot/hpbase/DAPT_cov/Landscaper/plot/{p}', p=PLOTFILES)
    output:
        'plot/hpbase/cont_DAPT_cov/energy_diff.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_cont_DAPT_cov.txt'
    log:
        'logs/featureplot_energy_cont_DAPT_cov.log'
    shell:
        'src/featureplot_energy_cont_DAPT_cov.sh >& {log}'

rule plot_cellular_density:
    input:
        'output/hpbase/{sample2}/seurat_annotated.RData'
    output:
        'plot/hpbase/{sample2}/cellular_density.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_cellular_density_{sample2}.txt'
    log:
        'logs/plot_cellular_density_{sample2}.log'
    shell:
        'src/plot_cellular_density.sh {input} {output} >& {log}'

rule dimplot_basin:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/{sample}/basin.png',
        'plot/hpbase/{sample}/basin_splitby.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
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
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_states_{sample}.txt'
    log:
        'logs/dimplot_states_{sample}.log'
    shell:
        'src/dimplot_states_{wildcards.sample}.sh {output} >& {log}'

rule plot_landscape:
    input:
        expand('plot/hpbase/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/hpbase/cont/landscape_rescaled.png',
        'plot/hpbase/DAPT/landscape_rescaled.png',
        'plot/hpbase/cont_cov/landscape_rescaled.png',
        'plot/hpbase/DAPT_cov/landscape_rescaled.png'
    container:
        'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_landscape.txt'
    log:
        'logs/plot_landscape.log'
    shell:
        'src/plot_landscape.sh {output} >& {log}'

rule plot_stable_states:
    input:
        'plot/hpbase/cont/Landscaper/Basin.tsv',
        'plot/hpbase/DAPT/Landscaper/Basin.tsv',
        'plot/hpbase/cont_cov/Landscaper/Basin.tsv',
        'plot/hpbase/DAPT_cov/Landscaper/Basin.tsv'
    output:
        'plot/hpbase/stable_states.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_stable_states.txt'
    log:
        'logs/plot_stable_states.log'
    shell:
        'src/plot_stable_states.sh {output} >& {log}'

