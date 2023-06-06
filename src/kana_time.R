source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]
outfile4 <- commandArgs(trailingOnly=TRUE)[5]
outfile5 <- commandArgs(trailingOnly=TRUE)[6]

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("-24h", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile1)

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("-36h", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile2)

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("-48h", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile3)

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("-72h", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile4)

# Loading
sce <- readRDS(infile)
# Stratify
sce <- sce[, grep("-96h", colData(sce)$sample)]
# Save
saveRDS(sce, file=outfile5)
