source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile1)
bindata <- read.table(infile2, header=FALSE)

# Stratification
idx_cont <- grep("cont-", seurat.integrated@meta.data$sample)
idx_dapt <- grep("DAPT-", seurat.integrated@meta.data$sample)
bindata_cont <- bindata[idx_cont, ]
bindata_dapt <- bindata[idx_dapt, ]

# Save
write.table(bindata_cont, outfile1, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(bindata_dapt, outfile2, row.names=FALSE, col.names=FALSE, quote=FALSE)
