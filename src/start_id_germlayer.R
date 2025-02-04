source("src/Functions3.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
# infile = "output/hpbase/integrated/seurat_annotated.RData"

# Loading
load(infile)

# Preprocess
mat <- seurat.integrated[["SCT"]]@data
target <- grep("24h", seurat.integrated@meta.data$sample)
mat <- mat[, target]
germlayer <- seurat.integrated@meta.data$germlayer[target]

mat@x <- log2(mat@x + 1)
mat <- t(mat)

# Filtering
start_id <- .nearest_to_centroid(mat, germlayer)

# Save
save(start_id, file=outfile)
