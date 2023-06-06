source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("cont-", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile1)

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("DAPT-", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile2)
