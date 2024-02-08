source("src/Functions.R")

# Parameter
set.seed(1234)
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]
outfile4 <- commandArgs(trailingOnly=TRUE)[5]
outfile5 <- commandArgs(trailingOnly=TRUE)[6]
outfile6 <- commandArgs(trailingOnly=TRUE)[7]
outfile7 <- commandArgs(trailingOnly=TRUE)[8]

# Loading
load(infile)

# Rename
tmp <- seurat.integrated

# Stratify
target.cont <- grep("cont-", tmp@meta.data$sample)
target.dapt <- grep("DAPT-", tmp@meta.data$sample)
target.24h <- grep("24h", tmp@meta.data$sample)
target.36h <- grep("36h", tmp@meta.data$sample)
target.48h <- grep("48h", tmp@meta.data$sample)
target.72h <- grep("72h", tmp@meta.data$sample)
target.96h <- grep("96h", tmp@meta.data$sample)

# Save
seurat.integrated <- tmp[, target.cont]
save(seurat.integrated, file=outfile1)

# Save
seurat.integrated <- tmp[, target.dapt]
save(seurat.integrated, file=outfile2)

# Save
seurat.integrated <- tmp[, target.24h]
save(seurat.integrated, file=outfile3)

# Save
seurat.integrated <- tmp[, target.36h]
save(seurat.integrated, file=outfile4)

# Save
seurat.integrated <- tmp[, target.48h]
save(seurat.integrated, file=outfile5)

# Save
seurat.integrated <- tmp[, target.72h]
save(seurat.integrated, file=outfile6)

# Save
seurat.integrated <- tmp[, target.96h]
save(seurat.integrated, file=outfile7)
