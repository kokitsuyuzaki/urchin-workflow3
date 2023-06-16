# DAG graph
snakemake -s workflow/download.smk --rulegraph | dot -Tpng > plot/download.png
snakemake -s workflow/gtf.smk --rulegraph | dot -Tpng > plot/gtf.png
snakemake -s workflow/preprocess.smk --rulegraph | dot -Tpng > plot/preprocess.png
snakemake -s workflow/cellranger.smk --rulegraph | dot -Tpng > plot/cellranger.png
snakemake -s workflow/seurat.smk --rulegraph | dot -Tpng > plot/seurat.png
snakemake -s workflow/doublet.smk --rulegraph | dot -Tpng > plot/doublet.png
snakemake -s workflow/trajectory.smk --rulegraph | dot -Tpng > plot/trajectory.png
snakemake -s workflow/velocity.smk --rulegraph | dot -Tpng > plot/velocity.png
snakemake -s workflow/plot.smk --rulegraph | dot -Tpng > plot/plot.png
snakemake -s workflow/report.smk --rulegraph | dot -Tpng > plot/report.png
snakemake -s workflow/celltyping.smk --rulegraph | dot -Tpng > plot/celltyping.png
snakemake -s workflow/template_matching.smk --rulegraph | dot -Tpng > plot/template_matching.png

