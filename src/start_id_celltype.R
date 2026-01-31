source("src/Functions3.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
# infile = 'output/hpbase/DAPT/seurat_annotated_landscaper.RData'

# Loading
load(infile)

# Preprocess
mat <- seurat.integrated[["SCT"]]@data
target <- grep("36h", seurat.integrated@meta.data$sample)
mat <- mat[, target]
celltype <- seurat.integrated@meta.data$celltype[target]

mat@x <- log2(mat@x + 1)
mat <- t(mat)

# Filtering
start_id <- .nearest_to_centroid(mat, celltype)

# Save
save(start_id, file=outfile)
