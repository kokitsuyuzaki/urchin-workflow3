source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly = TRUE)
infile <- args[1]
outfile <- args[2]

# Loading
load(infile)

# Pre-processing
df <- as.data.frame(seurat.integrated@reductions$umap@cell.embeddings)
colnames(df) <- c("UMAP_1", "UMAP_2")

# Plot
g <- ggplot(df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_density_2d_filled() + theme(legend.position = 'none')
ggsave(file=outfile, g, dpi=200, width=6, height=6)
