import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']
MODES = ['deterministic', 'stochastic', 'dynamical']
GERM_LAYER = ['ectoderm', 'mesoderm', 'endoderm']
rule all:
    input:
        expand('output/hpbase/integrated/velocyto/integrated_{mode}.h5ad',
            mode=MODES),
        expand('output/hpbase/cont/velocyto/cont_{mode}.h5ad',
            mode=MODES),
        expand('output/hpbase/DAPT/velocyto/DAPT_{mode}.h5ad',
            mode=MODES),
        expand('output/hpbase/cont/velocyto/cont_{gl}_{mode}.h5ad',
            gl=GERM_LAYER, mode=MODES),
        expand('output/hpbase/DAPT/velocyto/DAPT_{gl}_{mode}.h5ad',
            gl=GERM_LAYER, mode=MODES),
        expand('output/hpbase/integrated/velocyto/{gl}_{mode}.h5ad',
            gl=GERM_LAYER, mode=MODES)

#################################
# Generate Loom files
#################################
rule velocyto:
    input:
        'output/hpbase/{sample}/outs/web_summary.html',
        'data/hpbase/HpulGenome_v1_geneid.gtf',
        'output/echinobase/{sample}/outs/web_summary.html',
        'data/echinobase/sp5_0_GCF_geneid.gtf'
    output:
        'output/hpbase/{sample}/velocyto/{sample}.loom'
    wildcard_constraints:
        sample='|'.join([re.escape(x) for x in SAMPLES])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/velocyto_{sample}.txt'
    log:
        'logs/velocyto_{sample}.log'
    shell:
        'src/velocyto.sh {wildcards.sample} >& {log}'

#################################
# Aggregate Loom files
#################################
rule aggr_loom:
    input:
        expand('output/hpbase/{sample}/velocyto/{sample}.loom',
            sample=SAMPLES)
    output:
        'output/hpbase/aggr/velocyto/aggr.loom'
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/aggr_loom.txt'
    log:
        'logs/aggr_loom.log'
    shell:
        'src/aggr_loom.sh {output} >& {log}'

# rule aggr_loom_cont:
#     input:
#         expand('output/hpbase/{sample}/velocyto/{sample}.loom',
#             sample=SAMPLES)
#     output:
#         'output/hpbase/aggr_cont/velocyto/aggr.loom'
#     container:
#         'docker://koki/velocyto:20221005'
#     resources:
#         mem_mb=1000000
#     benchmark:
#         'benchmarks/aggr_loom_cont.txt'
#     log:
#         'logs/aggr_loom_cont.log'
#     shell:
#         'src/aggr_loom_cont.sh {output} >& {log}'

# rule aggr_loom_DAPT:
#     input:
#         expand('output/hpbase/{sample}/velocyto/{sample}.loom',
#             sample=SAMPLES)
#     output:
#         'output/hpbase/aggr_DAPT/velocyto/aggr.loom'
#     container:
#         'docker://koki/velocyto:20221005'
#     resources:
#         mem_mb=1000000
#     benchmark:
#         'benchmarks/aggr_loom_DAPT.txt'
#     log:
#         'logs/aggr_loom_DAPT.log'
#     shell:
#         'src/aggr_loom_DAPT.sh {output} >& {log}'

#################################
# Seurat => AnnData
#################################
rule seurat2anndata_integrated:
    input:
        'output/hpbase/integrated/seurat_annotated.RData'
    output:
        'output/hpbase/integrated/seurat_annotated.h5ad'
    container:
        'docker://koki/velocytor:20221015'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/seurat2anndata_integrated.txt'
    log:
        'logs/seurat2anndata_integrated.log'
    shell:
        'src/seurat2anndata_integrated.sh {input} {output} >& {log}'

#################################
# Calculte RNA Velocity
#################################
rule scvelo_integrated:
    input:
        'output/hpbase/aggr/velocyto/aggr.loom',
        'output/hpbase/integrated/seurat_annotated.h5ad'
    output:
        'output/hpbase/integrated/velocyto/integrated_{mode}.h5ad'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scvelo_integrated_{mode}.txt'
    log:
        'logs/scvelo_integrated_{mode}.log'
    shell:
        'src/scvelo.sh {wildcards.mode} {input} {output} >& {log}'

rule scvelo_germlayer:
    input:
        'output/hpbase/aggr/velocyto/aggr.loom',
        'output/hpbase/integrated/seurat_annotated.h5ad'
    output:
        'output/hpbase/integrated/velocyto/{gl}_{mode}.h5ad'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES]),
        gl='|'.join([re.escape(x) for x in GERM_LAYER])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scvelo_{gl}_{mode}.txt'
    log:
        'logs/scvelo_{gl}_{mode}.log'
    shell:
        'src/scvelo_germlayer.sh {wildcards.gl} {wildcards.mode} {input} {output} >& {log}'

rule scvelo_cont:
    input:
        'output/hpbase/aggr/velocyto/aggr.loom',
        'output/hpbase/integrated/seurat_annotated.h5ad'
    output:
        'output/hpbase/cont/velocyto/cont_{mode}.h5ad'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES]),
        gl='|'.join([re.escape(x) for x in GERM_LAYER])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scvelo_cont_{mode}.txt'
    log:
        'logs/scvelo_cont_{mode}.log'
    shell:
        'src/scvelo_cont.sh {wildcards.mode} {input} {output} >& {log}'

rule scvelo_DAPT:
    input:
        'output/hpbase/aggr/velocyto/aggr.loom',
        'output/hpbase/integrated/seurat_annotated.h5ad'
    output:
        'output/hpbase/DAPT/velocyto/DAPT_{mode}.h5ad'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES]),
        gl='|'.join([re.escape(x) for x in GERM_LAYER])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scvelo_DAPT_{mode}.txt'
    log:
        'logs/scvelo_DAPT_{mode}.log'
    shell:
        'src/scvelo_DAPT.sh {wildcards.mode} {input} {output} >& {log}'


rule scvelo_cont_germlayer:
    input:
        'output/hpbase/aggr/velocyto/aggr.loom',
        'output/hpbase/integrated/seurat_annotated.h5ad'
    output:
        'output/hpbase/cont/velocyto/cont_{gl}_{mode}.h5ad'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES]),
        gl='|'.join([re.escape(x) for x in GERM_LAYER])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scvelo_cont_{gl}_{mode}.txt'
    log:
        'logs/scvelo_cont_{gl}_{mode}.log'
    shell:
        'src/scvelo_germlayer.sh {wildcards.gl} {wildcards.mode} {input} {output} >& {log}'

rule scvelo_DAPT_germlayer:
    input:
        'output/hpbase/aggr/velocyto/aggr.loom',
        'output/hpbase/integrated/seurat_annotated.h5ad'
    output:
        'output/hpbase/DAPT/velocyto/DAPT_{gl}_{mode}.h5ad'
    wildcard_constraints:
        mode='|'.join([re.escape(x) for x in MODES]),
        gl='|'.join([re.escape(x) for x in GERM_LAYER])
    container:
        'docker://koki/velocyto:20221005'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/scvelo_DAPT_{gl}_{mode}.txt'
    log:
        'logs/scvelo_DAPT_{gl}_{mode}.log'
    shell:
        'src/scvelo_germlayer.sh {wildcards.gl} {wildcards.mode} {input} {output} >& {log}'
