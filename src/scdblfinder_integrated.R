source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Predict Doublets
sce <- as.SingleCellExperiment(seurat.integrated)
dbl.dens <- computeDoubletDensity(sce,  d=ncol(reducedDim(sce)))

# Save
save(dbl.dens, file=outfile)