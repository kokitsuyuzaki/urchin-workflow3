source("src/Functions2.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Conversion
outfile2 <- gsub(".h5ad", ".h5Seurat", outfile)
SaveH5Seurat(seurat.obj, filename = outfile2)
Convert(outfile2, dest = "h5ad")
