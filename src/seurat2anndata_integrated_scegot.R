source("src/Functions2.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# 1. sample列から'cont-'と'DAPT-'を削除
seurat.integrated$cluster_day <- gsub("^cont-", "", seurat.integrated$sample)
seurat.integrated$cluster_day <- gsub("^DAPT-", "", seurat.integrated$cluster_day)

# data slot を counts に置き換え
DefaultAssay(seurat.integrated) <- "RNA"
seurat.integrated[["RNA"]]@data <- seurat.integrated[["RNA"]]@counts

# Conversion
outfile2 <- gsub(".h5ad", ".h5Seurat", outfile)
unlink(outfile2)
SaveH5Seurat(seurat.integrated, filename = outfile2)
Convert(outfile2, dest = "h5ad")
