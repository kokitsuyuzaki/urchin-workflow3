import pandas as pd
from snakemake.utils import min_version

#################################
# Setting
#################################
min_version("6.5.3")

SAMPLES = ['cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h']

rule all:
    input:
        'plot/hpbase/integrated/proportion_plot_supervised.png',
        'plot/hpbase/integrated/proportion_plot_supervised_neurons.png',
        'plot/hpbase/integrated/proportion_plot_guidedpca.png',
        'plot/hpbase/integrated/proportion_plot_guidedpca_neurons.png',
        expand('plot/hpbase/{sample}/proportion_plot_guidedpca.png',
            sample=SAMPLES)

#################################
# Preprocessing
#################################
rule signalnoise_preprocess:
    input:
        'output/hpbase/integrated/seurat_annotated.RData',
        'output/hpbase/integrated/scdblfinder.RData',
        'data/annotation.RData'
    output:
        'output/hpbase/integrated/signalnoise.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_preprocess_hpbase_integrated.txt'
    log:
        'logs/signalnoise_preprocess_hpbase_integrated.log'
    shell:
        'src/signalnoise_preprocess.sh {input} {output} >& {log}'

# Stratification
rule signalnoise_label_integrated_neurons:
    input:
        'output/hpbase/integrated/signalnoise.RData'
    output:
        'output/hpbase/integrated/signalnoise_neurons.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_label_integrated_neurons.txt'
    log:
        'logs/signalnoise_label_integrated_neurons.log'
    shell:
        'src/label_integrated_neurons.sh {input} {output} >& {log}'

rule signalnoise_label_sample:
    input:
        'output/hpbase/integrated/signalnoise.RData'
    output:
        'output/hpbase/{sample}/signalnoise.RData'
    container:
        'docker://koki/urchin_workflow_seurat:20251014'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_label_integrated_{sample}.txt'
    log:
        'logs/signalnoise_label_integrated_{sample}.log'
    shell:
        'src/label_sample.sh {wildcards.sample} {input} {output} >& {log}'

#################################
# Supervised learning
#################################
rule signalnoise_supervised_learning:
    input:
        'output/hpbase/integrated/signalnoise.RData'
    output:
        'output/hpbase/integrated/signalnoise_supervised.RData'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_supervised_learning.txt'
    log:
        'logs/signalnoise_supervised_learning.log'
    shell:
        'src/signalnoise_supervised_learning.sh {input} {output} >& {log}'

rule signalnoise_supervised_learning_neurons:
    input:
        'output/hpbase/integrated/signalnoise_neurons.RData'
    output:
        'output/hpbase/integrated/signalnoise_supervised_neurons.RData'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_supervised_learning_neurons.txt'
    log:
        'logs/signalnoise_supervised_learning_neurons.log'
    shell:
        'src/signalnoise_supervised_learning.sh {input} {output} >& {log}'

#################################
# Guided PCA
#################################
rule signalnoise_guidedPCA:
    input:
        'output/hpbase/integrated/signalnoise.RData'
    output:
        'output/hpbase/integrated/signalnoise_guidedpca.RData'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_guidedPCA.txt'
    log:
        'logs/signalnoise_guidedPCA.log'
    shell:
        'src/signalnoise_guidedPCA.sh {input} {output} >& {log}'

rule signalnoise_guidedPCA_neurons:
    input:
        'output/hpbase/integrated/signalnoise_neurons.RData'
    output:
        'output/hpbase/integrated/signalnoise_guidedpca_neurons.RData'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_guidedPCA_neurons.txt'
    log:
        'logs/signalnoise_guidedPCA_neurons.log'
    shell:
        'src/signalnoise_guidedPCA.sh {input} {output} >& {log}'

rule signalnoise_guidedPCA_sample:
    input:
        'output/hpbase/{sample}/signalnoise.RData'
    output:
        'output/hpbase/{sample}/signalnoise_guidedpca.RData'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/signalnoise_guidedPCA_{sample}.txt'
    log:
        'logs/signalnoise_guidedPCA_{sample}.log'
    shell:
        'src/signalnoise_guidedPCA.sh {input} {output} >& {log}'

#################################
# Visualization
#################################
rule proportion_plot_signalnoise_supervised_learning:
    input:
        'output/hpbase/integrated/signalnoise_supervised.RData'
    output:
        'plot/hpbase/integrated/proportion_plot_supervised.png'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/proportion_plot_signalnoise_supervised_learning.txt'
    log:
        'logs/proportion_plot_signalnoise_supervised_learning.log'
    shell:
        'src/proportion_plot_signalnoise_supervised_learning.sh {input} {output} >& {log}'

rule proportion_plot_signalnoise_supervised_learning_neurons:
    input:
        'output/hpbase/integrated/signalnoise_supervised_neurons.RData'
    output:
        'plot/hpbase/integrated/proportion_plot_supervised_neurons.png'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/proportion_plot_signalnoise_supervised_learning_neurons.txt'
    log:
        'logs/proportion_plot_signalnoise_supervised_learning_neurons.log'
    shell:
        'src/proportion_plot_signalnoise_supervised_learning.sh {input} {output} >& {log}'

rule proportion_plot_signalnoise_guidedPCA:
    input:
        'output/hpbase/integrated/signalnoise_guidedpca.RData'
    output:
        'plot/hpbase/integrated/proportion_plot_guidedpca.png'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/proportion_plot_signalnoise_guidedPCA.txt'
    log:
        'logs/proportion_plot_signalnoise_guidedPCA.log'
    shell:
        'src/proportion_plot_signalnoise_guidedPCA.sh {input} {output} >& {log}'

rule proportion_plot_signalnoise_guidedPCA_neurons:
    input:
        'output/hpbase/integrated/signalnoise_guidedpca_neurons.RData'
    output:
        'plot/hpbase/integrated/proportion_plot_guidedpca_neurons.png'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/proportion_plot_signalnoise_guidedPCA_neurons.txt'
    log:
        'logs/proportion_plot_signalnoise_guidedPCA_neurons.log'
    shell:
        'src/proportion_plot_signalnoise_guidedPCA.sh {input} {output} >& {log}'

rule proportion_plot_signalnoise_guidedPCA_sample:
    input:
        'output/hpbase/{sample}/signalnoise_guidedpca.RData'
    output:
        'plot/hpbase/{sample}/proportion_plot_guidedpca.png'
    container:
        'docker://koki/urchin_workflow_signalnoise:20250821'
    resources:
        mem_mb=1000000
    benchmark:
        'benchmarks/proportion_plot_signalnoise_guidedPCA_{sample}.txt'
    log:
        'logs/proportion_plot_signalnoise_guidedPCA_{sample}.log'
    shell:
        'src/proportion_plot_signalnoise_guidedPCA_sample.sh {input} {output} >& {log}'
