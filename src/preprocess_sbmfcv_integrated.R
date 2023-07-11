source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Preprocess
out <- as.matrix(t(seurat.integrated@assays$RNA@counts[seurat.integrated@assays$integrated@var.features, ]))

# Output
write.table(out, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)