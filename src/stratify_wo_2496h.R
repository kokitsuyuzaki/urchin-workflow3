source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Stratification
idx <- which(seurat.integrated@meta.data$sample %in% c("cont-36h", "cont-48h", "cont-72h", "DAPT-36h", "DAPT-48h", "DAPT-72h"))
seurat.obj <- seurat.integrated[, idx]

# Save
save(seurat.obj, file=outfile)
