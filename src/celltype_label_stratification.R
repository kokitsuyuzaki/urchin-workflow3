source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# Stratify
target.cont <- grep("cont-", seurat.integrated@meta.data$sample)
target.dapt <- grep("DAPT-", seurat.integrated@meta.data$sample)
seurat.cont <- seurat.integrated[, target.cont]
seurat.dapt <- seurat.integrated[, target.dapt]

# Save
save(seurat.cont, file=outfile1)
save(seurat.dapt, file=outfile2)
