source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Setting for kana
sce <- as.SingleCellExperiment(seurat.obj)
counts(sce) <- NULL
altExp(sce) <- NULL

# Save
saveRDS(sce, file=outfile)
