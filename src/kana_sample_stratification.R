source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
sample <- commandArgs(trailingOnly=TRUE)[3]

# Loading
sce <- readRDS(infile)

# Stratification
sce <- sce[, which(colData(sce)$sample == sample)]

# Save
saveRDS(sce, file=outfile)
