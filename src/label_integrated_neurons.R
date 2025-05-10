source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Sub-sampling
seurat.integrated <- subset(seurat.integrated, subset = celltype == "Neurons")
seurat.integrated <- subset(seurat.integrated, subset = sample %in% c("cont-72h", "cont-96h", "DAPT-72h", "DAPT-96h"))

# Re-clustering
seurat.integrated <- FindNeighbors(seurat.integrated, dims=1:30)
seurat.integrated <- FindClusters(seurat.integrated)


# Save
save(seurat.integrated, file=outfile)
