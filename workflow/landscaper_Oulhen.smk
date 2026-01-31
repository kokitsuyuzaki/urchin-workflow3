import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ["integrated", "cont", "DAPT"]
SAMPLES2 = ["integrated", "cont_stratified", "DAPT_stratified"]
PLOTFILES = ['ratio_group.png', 'Allstates.png', 'Freq_Prob_Energy.png', 'h.png', 'J.png', 'Basin.png', 'StatusNetwork_Subgraph.png', 'StatusNetwork_Subgraph_legend.png', 'StatusNetwork_Energy.png', 'StatusNetwork_Energy_legend.png', 'StatusNetwork_Ratio.png', 'StatusNetwork_Ratio_legend.png', 'StatusNetwork_State.png', 'StatusNetwork_State_legend.png', 'Landscape.png', 'discon_graph_1.png', 'discon_graph_2.png']

rule all:
    input:
        'plot/echinobase/Oulhen/integrated/pca_landscaper.png',
        expand('plot/echinobase/Oulhen/{sample}/dimplot_celltype_landscaper.png',
            sample=SAMPLES),
        'plot/echinobase/Oulhen/integrated/featureplot_ncount_rna_landscaper.png',
        # # expand('plot/echinobase/Oulhen/{sample}/binpca/deg/FINISH',
        # #     sample=SAMPLES),        
        # expand('plot/echinobase/Oulhen/{sample}/bindata/FINISH',
        #     sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/sce.RData',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/energy.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/energy_hex.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/energy_rescaled.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/energy_rescaled_hex.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/energy_contour.png',
            sample=SAMPLES),
        'plot/echinobase/Oulhen/cont_DAPT/energy_diff.png',
        expand('plot/echinobase/Oulhen/{sample}/cellular_density.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/basin.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/states.png',
            sample=SAMPLES),
        # 'plot/echinobase/Oulhen/cont/landscape_rescaled.png',
        # 'plot/echinobase/Oulhen/DAPT/landscape_rescaled.png',
        # 'plot/echinobase/Oulhen/cont_cov/landscape_rescaled.png',
        # 'plot/echinobase/Oulhen/DAPT_cov/landscape_rescaled.png',
        # Basin Piechart
        expand('plot/echinobase/Oulhen/{sample}/basin_piechart.png',
            sample=SAMPLES),
        # Transition Probability
        expand('plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Aboral_ectoderm.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Oral_ectoderm.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Ciliary_band.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Neurons.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_glauber_rw_Aboral_ectoderm.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_glauber_rw_Oral_ectoderm.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_glauber_rw_Ciliary_band.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_glauber_rw_Neurons.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/Landscaper/plot/discon_graph_1_new.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/Landscaper/plot/discon_graph_2_new.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_metropolis_emded.png',
            sample=SAMPLES),
        expand('plot/echinobase/Oulhen/{sample}/P_glauber_emded.png',
            sample=SAMPLES),
        # # Coarse-grained Vector Fields
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_germlayer.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_cluster.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_sample.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_celltype.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_state.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_energy.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_germlayer.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_cluster.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_sample.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_celltype.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_state.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_energy.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_germlayer_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_cluster_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_sample_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_celltype_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_state_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_metropolis_cg_energy_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_germlayer_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_cluster_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_sample_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_celltype_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_state_streamline.png',
        #     sample=SAMPLES),
        # expand('plot/echinobase/Oulhen/{sample}/P_glauber_cg_energy_streamline.png',
        #     sample=SAMPLES),

#######################################
# Dimension Reduction & Binarization
#######################################
rule preprocess_landscaper_Oulhen:
    input:
        'data/Oulhen/Combined_Dim20_res0.5.rds'
    output:
        'output/echinobase/Oulhen/integrated/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/integrated/seurat.tsv',
        'output/echinobase/Oulhen/integrated/group.tsv',
        'output/echinobase/Oulhen/cont/group.tsv',
        'output/echinobase/Oulhen/DAPT/group.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/preprocess_landscaper_Oulhen.txt'
    log:
        'logs/preprocess_landscaper_Oulhen.log'
    shell:
        'src/preprocess_landscaper_Oulhen.sh >& {log}'

rule binpca:
    input:
        'output/echinobase/Oulhen/integrated/seurat.tsv'
    output:
        'output/echinobase/Oulhen/integrated/binpca/BIN_DATA.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/binpca.txt'
    log:
        'logs/binpca.log'
    shell:
        'src/binpca.sh {input} {output} >& {log}'

rule pca_landscaper_Oulhen:
    input:
        'output/echinobase/Oulhen/integrated/seurat_annotated_landscaper.RData',
    output:
        'plot/echinobase/Oulhen/integrated/pca_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/pca_landscaper_Oulhen.txt'
    log:
        'logs/pca_landscaper_Oulhen.log'
    shell:
        'src/pca_landscaper.sh {input} {output} >& {log}'

rule stratify_Oulhen:
    input:
        'output/echinobase/Oulhen/integrated/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/integrated/binpca/BIN_DATA.tsv'
    output:
        'output/echinobase/Oulhen/cont/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/DAPT/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/cont/binpca/BIN_DATA.tsv',
        'output/echinobase/Oulhen/DAPT/binpca/BIN_DATA.tsv'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/stratify_Oulhen.txt'
    log:
        'logs/stratify_Oulhen.log'
    shell:
        'src/stratify_Oulhen.sh {input} {output} >& {log}'

rule dimplot_celltype_landscaper_Oulhen:
    input:
        'output/echinobase/Oulhen/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/echinobase/Oulhen/{sample}/dimplot_celltype_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/dimplot_celltype_landscaper_Oulhen_{sample}.txt'
    log:
        'logs/dimplot_celltype_landscaper_Oulhen_{sample}.log'
    shell:
        'src/dimplot_celltype_landscaper_Oulhen.sh {input} {output} >& {log}'

rule featureplot_ncount_rna_landscaper_Oulhen:
    input:
        'output/echinobase/Oulhen/integrated/seurat_annotated_landscaper.RData'
    output:
        'plot/echinobase/Oulhen/integrated/featureplot_ncount_rna_landscaper.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/featureplot_ncount_rna_landscaper_Oulhen.txt'
    log:
        'logs/featureplot_ncount_rna_landscaper_Oulhen.log'
    shell:
        'src/featureplot_ncount_rna_landscaper_Oulhen.sh {input} {output} >& {log}'

rule dimplot_bindata:
    input:
        'output/echinobase/Oulhen/{sample}/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/{sample}/binpca/BIN_DATA.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/bindata/FINISH'
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
        'src/dimplot_bindata_{wildcards.sample}.sh {input} {output} >& {log}'

# rule deg_bindata:
#     input:
#         'output/echinobase/{sample}/seurat_annotated.RData',
#         'output/echinobase/{sample}/binpca/BIN_DATA.tsv'
#     output:
#         'output/echinobase/{sample}/binpca/deg_bindata.xlsx'
#     container:
#         'docker://koki/urchin_workflow_seurat:20251014'
#     resources:
#         mem_mb=10000000
#     benchmark:
#         'benchmarks/deg_bindata_{sample}.txt'
#     log:
#         'logs/deg_bindata_{sample}.log'
#     shell:
#         'src/deg_bindata.sh {input} {output} >& {log}'

# rule featureplot_deg_bindata:
#     input:
#         'output/echinobase/{sample}/seurat_annotated.RData',
#         'output/echinobase/{sample}/binpca/deg_bindata.xlsx'
#     output:
#         'plot/echinobase/{sample}/binpca/deg/FINISH'
#     container:
#         'docker://koki/urchin_workflow_seurat:20251014'
#     resources:
#         mem_gb=1000
#     benchmark:
#         'benchmarks/featureplot_deg_bindata_{sample}.txt'
#     log:
#         'logs/featureplot_deg_bindata_{sample}.log'
#     shell:
#         'src/featureplot_deg_bindata_{wildcards.sample}.sh {input} {output} >& {log}'

rule landscaper_integrated:
    input:
        'output/echinobase/Oulhen/integrated/binpca/BIN_DATA.tsv',
        'output/echinobase/Oulhen/integrated/group.tsv'
    output:
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/ratio_group.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/Allstates.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/h.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/J.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/Basin.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_State.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/Landscape.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/discon_graph_1.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/plot/discon_graph_2.png',
        'plot/echinobase/Oulhen/integrated/Landscaper/Coordinate.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/Allstates.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/Basin.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/major_group.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/igraph.RData',
        'plot/echinobase/Oulhen/integrated/Landscaper/EnergyBarrier.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/Allstates_major_group.tsv'
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
        'output/echinobase/Oulhen/cont/binpca/BIN_DATA.tsv',
        'output/echinobase/Oulhen/cont/group.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/Coordinate.tsv'
    output:
        'plot/echinobase/Oulhen/cont/Landscaper/plot/ratio_group.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/Allstates.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/h.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/J.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/Basin.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_State.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/Landscape.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/discon_graph_1.png',
        'plot/echinobase/Oulhen/cont/Landscaper/plot/discon_graph_2.png',
        'plot/echinobase/Oulhen/cont/Landscaper/Basin.tsv',
        'plot/echinobase/Oulhen/cont/Landscaper/major_group.tsv',
        'plot/echinobase/Oulhen/cont/Landscaper/Allstates.tsv',
        'plot/echinobase/Oulhen/cont/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/cont/Landscaper/igraph.RData',
        'plot/echinobase/Oulhen/cont/Landscaper/EnergyBarrier.tsv',
        'plot/echinobase/Oulhen/cont/Landscaper/Allstates_major_group.tsv'
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
        'output/echinobase/Oulhen/DAPT/binpca/BIN_DATA.tsv',
        'output/echinobase/Oulhen/DAPT/group.tsv',
        'plot/echinobase/Oulhen/integrated/Landscaper/Coordinate.tsv'
    output:
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/ratio_group.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/Allstates.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/Freq_Prob_Energy.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/h.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/J.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/Basin.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_Subgraph.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_Subgraph_legend.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_Energy.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_Energy_legend.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_Ratio.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_Ratio_legend.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_State.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/StatusNetwork_State_legend.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/Landscape.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/discon_graph_1.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/plot/discon_graph_2.png',
        'plot/echinobase/Oulhen/DAPT/Landscaper/Basin.tsv',
        'plot/echinobase/Oulhen/DAPT/Landscaper/major_group.tsv',
        'plot/echinobase/Oulhen/DAPT/Landscaper/Allstates.tsv',
        'plot/echinobase/Oulhen/DAPT/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/DAPT/Landscaper/igraph.RData',
        'plot/echinobase/Oulhen/DAPT/Landscaper/EnergyBarrier.tsv',
        'plot/echinobase/Oulhen/DAPT/Landscaper/Allstates_major_group.tsv'
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

rule transition_probability_Oulhen:
    input:
        'plot/echinobase/Oulhen/{sample}/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Allstates.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/P_metropolis.tsv',
        'plot/echinobase/Oulhen/{sample}/P_glauber.tsv'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/transition_probability_Oulhen_{sample}.txt'
    log:
        'logs/transition_probability_Oulhen_{sample}.log'
    shell:
        'src/transition_probability_Oulhen.sh {input} {output} >& {log}'

rule coarse_grained_vector_fields:
    input:
        'output/echinobase/Oulhen/integrated/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/{sample}/seurat_annotated_landscaper.RData',
        'output/echinobase/Oulhen/{sample}/binpca/BIN_DATA.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Allstates.tsv',
        'plot/echinobase/Oulhen/{sample}/P_metropolis.tsv',
        'plot/echinobase/Oulhen/{sample}/P_glauber.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/P_metropolis_cg.RData',
        'plot/echinobase/Oulhen/{sample}/P_glauber_cg.RData'
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

rule random_walk:
    input:
        'plot/echinobase/Oulhen/{sample}/P_metropolis.tsv',
        'plot/echinobase/Oulhen/{sample}/P_glauber.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/P_metropolis_rw.tsv',
        'plot/echinobase/Oulhen/{sample}/P_glauber_rw.tsv'
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

rule plot_random_walk_metropolis_Oulhen:
    input:
        'plot/echinobase/Oulhen/integrated/Landscaper/major_group.tsv',
        'plot/echinobase/Oulhen/{sample}/P_metropolis_rw.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Aboral_ectoderm.png',
        'plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Oral_ectoderm.png',
        'plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Ciliary_band.png',
        'plot/echinobase/Oulhen/{sample}/P_metropolis_rw_Neurons.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_metropolis_Oulhen_{sample}.txt'
    log:
        'logs/plot_random_walk_metropolis_Oulhen_{sample}.log'
    shell:
        'src/plot_random_walk_Oulhen.sh {input} {output} >& {log}'

rule plot_random_walk_glauber_Oulhen:
    input:
        'plot/echinobase/Oulhen/integrated/Landscaper/major_group.tsv',
        'plot/echinobase/Oulhen/{sample}/P_glauber_rw.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/P_glauber_rw_Aboral_ectoderm.png',
        'plot/echinobase/Oulhen/{sample}/P_glauber_rw_Oral_ectoderm.png',
        'plot/echinobase/Oulhen/{sample}/P_glauber_rw_Ciliary_band.png',
        'plot/echinobase/Oulhen/{sample}/P_glauber_rw_Neurons.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_random_walk_glauber_Oulhen_{sample}.txt'
    log:
        'logs/plot_random_walk_glauber_Oulhen_{sample}.log'
    shell:
        'src/plot_random_walk_Oulhen.sh {input} {output} >& {log}'

rule graph_embedding:
    input:
        'plot/echinobase/Oulhen/{sample}/Landscaper/Allstates.tsv',
        'plot/echinobase/Oulhen/{sample}/P_metropolis.tsv',
        'plot/echinobase/Oulhen/{sample}/P_glauber.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/P_metropolis_emded.RData',
        'plot/echinobase/Oulhen/{sample}/P_glauber_emded.RData'
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
        'plot/echinobase/Oulhen/{sample}/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Basin.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/igraph.RData',
        'plot/echinobase/Oulhen/{sample}/P_metropolis_emded.RData',
        'plot/echinobase/Oulhen/{sample}/P_glauber_emded.RData'
    output:
        'plot/echinobase/Oulhen/{sample}/P_metropolis_emded.png',
        'plot/echinobase/Oulhen/{sample}/P_glauber_emded.png'
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
        'plot/echinobase/Oulhen/{sample}/Landscaper/Allstates.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Basin.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/Landscaper/TippingState.tsv'
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
        'plot/echinobase/Oulhen/{sample}/Landscaper/E.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Basin.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/EnergyBarrier.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Allstates_major_group.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/TippingState.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/Landscaper/dendrogram_new.RData'
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
        'plot/echinobase/Oulhen/{sample}/Landscaper/dendrogram_new.RData'
    output:
        'plot/echinobase/Oulhen/{sample}/Landscaper/plot/discon_graph_1_new.png',
        'plot/echinobase/Oulhen/{sample}/Landscaper/plot/discon_graph_2_new.png'
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

# rule plot_h:
#     input:
#         'plot/echinobase/Oulhen/cont/Landscaper/plot/h.png',
#         'plot/echinobase/Oulhen/DAPT/Landscaper/plot/h.png',
#     output:
#         'plot/echinobase/Oulhen/integrated/h.png'
#     container:
#         'docker://koki/urchin_workflow_seurat:20251014'
#     resources:
#         mem_mb=10000000
#     benchmark:
#         'benchmarks/plot_h.txt'
#     log:
#         'logs/plot_h.log'
#     shell:
#         'src/plot_h.sh {output} >& {log}'

# rule plot_J:
#     input:
#         'plot/echinobase/Oulhen/cont/Landscaper/plot/J.png',
#         'plot/echinobase/Oulhen/DAPT/Landscaper/plot/J.png'
#     output:
#         'plot/echinobase/Oulhen/integrated/J.png'
#     container:
#         'docker://koki/urchin_workflow_seurat:20251014'
#     resources:
#         mem_mb=10000000
#     benchmark:
#         'benchmarks/plot_J.txt'
#     log:
#         'logs/plot_J.log'
#     shell:
#         'src/plot_J.sh {output} >& {log}'

rule featureplot_energy_Oulhen:
    input:
        expand('plot/echinobase/Oulhen/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/echinobase/Oulhen/{sample}/sce.RData',
        'plot/echinobase/Oulhen/{sample}/energy.png',
        'plot/echinobase/Oulhen/{sample}/energy_hex.png',
        'plot/echinobase/Oulhen/{sample}/energy_rescaled.png',
        'plot/echinobase/Oulhen/{sample}/energy_rescaled_hex.png',
        'plot/echinobase/Oulhen/{sample}/energy_contour.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_Oulhen_{sample}.txt'
    log:
        'logs/featureplot_energy_Oulhen_{sample}.log'
    shell:
        'src/featureplot_energy_Oulhen_{wildcards.sample}.sh {output} >& {log}'

rule featureplot_energy_Oulhen_cont_dapt:
    input:
        'plot/echinobase/Oulhen/integrated/sce.RData',
        expand('plot/echinobase/Oulhen/cont/Landscaper/plot/{p}', p=PLOTFILES),
        expand('plot/echinobase/Oulhen/DAPT/Landscaper/plot/{p}', p=PLOTFILES)
    output:
        'plot/echinobase/Oulhen/cont_DAPT/energy_diff.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/featureplot_energy_Oulhen_cont_DAPT.txt'
    log:
        'logs/featureplot_energy_Oulhen_cont_DAPT.log'
    shell:
        'src/featureplot_energy_Oulhen_cont_DAPT.sh >& {log}'

rule plot_cellular_density:
    input:
        'output/echinobase/Oulhen/{sample}/seurat_annotated_landscaper.RData'
    output:
        'plot/echinobase/Oulhen/{sample}/cellular_density.png'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/plot_cellular_density_{sample}.txt'
    log:
        'logs/plot_cellular_density_{sample}.log'
    shell:
        'src/plot_cellular_density.sh {input} {output} >& {log}'

rule dimplot_basin_Oulhen:
    input:
        expand('plot/echinobase/Oulhen/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/echinobase/Oulhen/{sample}/basin.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_basin_Oulhen_{sample}.txt'
    log:
        'logs/dimplot_basin_Oulhen_{sample}.log'
    shell:
        'src/dimplot_basin_Oulhen_{wildcards.sample}.sh {output} >& {log}'

rule piechart_basin_Oulhen:
    input:
        'output/echinobase/Oulhen/{sample}/seurat_annotated_landscaper.RData',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Allstates_major_group.tsv',
        'plot/echinobase/Oulhen/{sample}/Landscaper/BIN_DATA',
        'plot/echinobase/Oulhen/{sample}/Landscaper/Basin.tsv'
    output:
        'plot/echinobase/Oulhen/{sample}/basin_piechart.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/piechart_basin_Oulhen_{sample}.txt'
    log:
        'logs/piechart_basin_Oulhen_{sample}.log'
    shell:
        'src/piechart_basin_Oulhen.sh {input} {output} >& {log}'

rule dimplot_states_Oulhen:
    input:
        expand('plot/echinobase/Oulhen/{sample}/Landscaper/plot/{p}',
            sample=SAMPLES, p=PLOTFILES)
    output:
        'plot/echinobase/Oulhen/{sample}/states.png'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=10000000
    benchmark:
        'benchmarks/dimplot_states_Oulhen_{sample}.txt'
    log:
        'logs/dimplot_states_Oulhen_{sample}.log'
    shell:
        'src/dimplot_states_Oulhen_{wildcards.sample}.sh {output} >& {log}'

# rule plot_landscape:
#     input:
#         expand('plot/echinobase/Oulhen/{sample}/Landscaper/plot/{p}',
#             sample=SAMPLES, p=PLOTFILES)
#     output:
#         'plot/echinobase/Oulhen/cont/landscape_rescaled.png',
#         'plot/echinobase/Oulhen/DAPT/landscape_rescaled.png'
#     container:
#         'docker://ghcr.io/chiba-ai-med/landscaper:pr-26'
#     resources:
#         mem_mb=10000000
#     benchmark:
#         'benchmarks/plot_landscape.txt'
#     log:
#         'logs/plot_landscape.log'
#     shell:
#         'src/plot_landscape.sh {output} >& {log}'

# rule plot_stable_states:
#     input:
#         'plot/echinobase/Oulhen/cont/Landscaper/Basin.tsv',
#         'plot/echinobase/Oulhen/DAPT/Landscaper/Basin.tsv'
#     output:
#         'plot/echinobase/Oulhen/stable_states.png'
#     container:
#         'docker://koki/urchin_workflow_seurat:20251014'
#     resources:
#         mem_mb=10000000
#     benchmark:
#         'benchmarks/plot_stable_states.txt'
#     log:
#         'logs/plot_stable_states.log'
#     shell:
#         'src/plot_stable_states.sh {output} >& {log}'

