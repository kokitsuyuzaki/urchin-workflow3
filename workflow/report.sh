# HTML
mkdir -p report
snakemake -s workflow/download.smk --report report/download.html
snakemake -s workflow/gtf.smk --report report/gtf.html
snakemake -s workflow/preprocess.smk --report report/preprocess.html
snakemake -s workflow/cellranger.smk --report report/cellranger.html
snakemake -s workflow/seurat.smk --report report/seurat.html
snakemake -s workflow/doublet.smk --report report/doublet.html
snakemake -s workflow/trajectory.smk --report report/trajectory.html
snakemake -s workflow/velocity.smk --report report/velocity.html
snakemake -s workflow/plot.smk --report report/plot.html
snakemake -s workflow/report.smk --report report/report.html
snakemake -s workflow/celltyping.smk --report report/celltyping.html
snakemake -s workflow/template_matching.smk --report report/template_matching.html
snakemake -s workflow/landscaper.smk --report report/landscaper.html
snakemake -s workflow/signalnoise.smk --report report/signalnoise.html