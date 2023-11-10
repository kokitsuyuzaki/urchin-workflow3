source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Preprocess
## UMAP
# ## PCA（だめだったので消す）
# out <- seurat.integrated@reductions$pca@cell.embeddings[,1:20]
out <- RunUMAP(seurat.integrated, dims=1:50, n.components=20L)@reductions$umap@cell.embeddings

## Normalization
# Min-Max Normalization（だめだったので消す）
# out <- out - min(out)
out <- apply(out, 2, function(x){x - min(x)})

# Output
write.table(out, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)