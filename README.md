# urchin-workflow3
This workflow consists of 14 workflows as follows:

- **workflow/download.smk**: Data downloading

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/download.png?raw=true)

- **workflow/gtf.smk**: GTF file preprocessing

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/gtf.png?raw=true)

- **workflow/preprocess.smk**: Data preprocessing

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/preprocess.png?raw=true)

- **workflow/cellranger.smk**: `CellRanger` mkref and count

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/cellranger.png?raw=true)

- **workflow/seurat.smk**: `Seurat` against each dataset, data integration, and marker detection

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/seurat.png?raw=true)

- **workflow/doublet.smk**: Doublet detection by `scDblFinder`

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/doublet.png?raw=true)

- **workflow/trajectory.smk**: Trajectory inference by `Monocle3`

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/trajectory.png?raw=true)

- **workflow/plot.smk**: Plots against cluster labels (before celltyping)

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/plot.png?raw=true)

- **workflow/report.smk**: Interactive report by `scTGIF`

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/report.png?raw=true)

- **workflow/celltyping.smk**: Celltyping based on manual selection of cluster-specific markers

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/celltyping.png?raw=true)

- **workflow/velocity.smk**: RNA velocity analysis by `Velocyto` and `scVelo`.

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/velocity.png?raw=true)

- **workflow/template_matching.smk**: Template matching for Hp-Opn5L, Hp-Tph, and Hp-Delta

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/template_matching.png?raw=true)

- **workflow/landscaper.smk**: Energy Landscape Analysis by Landscaper

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/landscaper.png?raw=true)

- **workflow/signal_noise.smk**: Supervised Learning and Guided PCA

![](https://github.com/kokitsuyuzaki/urchin-workflow3/blob/main/plot/signal_noise.png?raw=true)

## Requirements
- Bash: GNU bash, version 4.2.46(1)-release (x86_64-redhat-linux-gnu)
- Snakemake: 6.5.3
- Singularity: 3.8.0

## How to reproduce this workflow
### In Local Machine

```
snakemake -s workflow/download.smk -j 4 --use-singularity
snakemake -s workflow/gtf.smk -j 4 --use-singularity
snakemake -s workflow/preprocess.smk -j 4 --use-singularity
snakemake -s workflow/cellranger.smk -j 4 --use-singularity
snakemake -s workflow/seurat.smk -j 4 --use-singularity
snakemake -s workflow/doublet.smk -j 4 --use-singularity
snakemake -s workflow/trajectory.smk -j 4 --use-singularity
snakemake -s workflow/velocity.smk -j 4 --use-singularity
snakemake -s workflow/plot.smk -j 4 --use-singularity
snakemake -s workflow/report.smk -j 4 --use-singularity
snakemake -s workflow/celltyping.smk -j 4 --use-singularity
snakemake -s workflow/template_matching.smk -j 4 --use-singularity
snakemake -s workflow/landscaper.smk -j 4 --use-singularity
snakemake -s workflow/signalnoise.smk -j 4 --use-singularity
```

### In Open Grid Engine

```
snakemake -s workflow/download.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/gtf.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/preprocess.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/cellranger.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/seurat.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/doublet.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/trajectory.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/velocity.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/plot.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/report.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/celltyping.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/template_matching.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/landscaper.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
snakemake -s workflow/signalnoise.smk -j 32 --cluster qsub --latency-wait 600 --use-singularity
```

### In Slurm

```
snakemake -s workflow/download.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/gtf.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/preprocess.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/cellranger.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/seurat.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/doublet.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/trajectory.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/velocity.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/plot.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/report.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/celltyping.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/template_matching.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/landscaper.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
snakemake -s workflow/signalnoise.smk -j 32 --cluster sbatch --latency-wait 600 --use-singularity
```

## License
Copyright (c) 2023 Koki Tsuyuzaki released under the [Artistic License 2.0](http://www.perlfoundation.org/artistic_license_2_0).

## Authors
- Koki Tsuyuzaki
- Shunsuke Yaguchi