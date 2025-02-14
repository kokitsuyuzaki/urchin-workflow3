source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile4 <- commandArgs(trailingOnly=TRUE)[4]
outfile5 <- commandArgs(trailingOnly=TRUE)[5]
outfile6 <- commandArgs(trailingOnly=TRUE)[6]

# Loading
all_states_cont <- unlist(read.delim('plot/hpbase/cont_cov/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/hpbase/cont/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_cont <- unlist(read.table('plot/hpbase/cont_cov/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/cont_stratified/seurat.RData')

# Sort
names(all_states_cont) <- NULL
names(bin_data_cont) <- NULL

target_cont <- sapply(bin_data_cont, function(x){
     which(all_states_cont == x)
})

energy_cont_sorted <- energy_cont[target_cont]
energy_cont_sorted <- as.matrix(energy_cont_sorted)
rownames(energy_cont_sorted) <- colnames(seurat.integrated)

# Assign Labels
seurat.integrated$energy <- energy_cont_sorted

# Seurat => SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = seurat.integrated@assays$RNA@counts))
reducedDims(sce) <- list(UMAP = seurat.integrated@reductions$umap@cell.embeddings)
logcounts(sce) <- log10(counts(sce) + 1)

# Assign Energy
colData(sce)$energy <- seurat.integrated$energy

# Setting Hex bin
sce <- make_hexbin(sce, nbins = 50, dimension_reduction = "UMAP")
save(sce, file=outfile1)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100))
ggsave(file=outfile2, g, dpi=200, width=6, height=6)

# Plot
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "")
ggsave(file=outfile3, g, dpi=200, width=6, height=6)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100), limits=c(-7.0, 8.5))
ggsave(file=outfile4, g, dpi=200, width=6, height=6)

# Plot
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "")
ggsave(file=outfile5, g, dpi=200, width=6, height=6)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", split.by="sample", ncol=5, pt.size=2)
png(file=outfile6, width=3000, height=600)
print(g)
dev.off()
