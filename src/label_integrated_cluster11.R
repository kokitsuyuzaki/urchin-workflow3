source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Sub-sampling
seurat.integrated <- subset(seurat.integrated, subset = seurat_clsters == "11")

# Re-clustering
seurat.integrated <- FindNeighbors(seurat.integrated, dims=1:30)
seurat.integrated <- FindClusters(seurat.integrated)

# Save
save(seurat.integrated, file=outfile)
