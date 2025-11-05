source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]
outfile4 <- commandArgs(trailingOnly=TRUE)[6]
outfile5 <- commandArgs(trailingOnly=TRUE)[7]
outfile6 <- commandArgs(trailingOnly=TRUE)[8]
outfile7 <- commandArgs(trailingOnly=TRUE)[9]
outfile8 <- commandArgs(trailingOnly=TRUE)[10]
outfile9 <- commandArgs(trailingOnly=TRUE)[11]
outfile10 <- commandArgs(trailingOnly=TRUE)[12]

# Loading
load(infile1)

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
seurat.integrated <- seurat.integrated[, intersect(target1, target2)]

seurat.integrated.original <- seurat.integrated
bindata <- read.table(infile2, header=FALSE)

# Stratification
idx_cont <- grep("cont-", seurat.integrated.original@meta.data$sample)
idx_dapt <- grep("DAPT-", seurat.integrated.original@meta.data$sample)

# Save（Seurat Object）
file.copy(infile1, outfile1)

seurat.integrated <- seurat.integrated.original[, idx_cont]
save(seurat.integrated, file=outfile2)
save(seurat.integrated, file=outfile3)

seurat.integrated <- seurat.integrated.original[, idx_dapt]
save(seurat.integrated, file=outfile4)
save(seurat.integrated, file=outfile5)

# Save（BIN Data）
file.copy(infile2, outfile6)

bindata_cont <- bindata[idx_cont, ]
bindata_dapt <- bindata[idx_dapt, ]

write.table(bindata_cont, outfile7, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bindata_cont, outfile8, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(bindata_dapt, outfile9, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bindata_dapt, outfile10, row.names=FALSE, col.names=FALSE, quote=FALSE)
