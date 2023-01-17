source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Log10
seurat.obj[["Log10_nCount_RNA"]] <- log10(seurat.obj@meta.data$nCount_RNA)

# Plot
png(file=outfile, width=600, height=600)
FeaturePlot(seurat.obj, features="Log10_nCount_RNA", pt.size=2, label.size=6)
dev.off()
