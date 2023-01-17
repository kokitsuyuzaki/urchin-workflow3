source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# Log10
seurat.integrated[["Log10_nCount_RNA"]] <- log10(seurat.integrated@meta.data$nCount_RNA)

# Plot
png(file=outfile1, width=600, height=600)
FeaturePlot(seurat.integrated, features="Log10_nCount_RNA",
    reduction = "umap", pt.size=2, label.size=6) + xlim(c(-15,15)) + ylim(c(-15,15))
dev.off()

# Plot
seuratList <- .stratifySeurat(seurat.integrated, group_names)
png(file=outfile2, width=2000, height=1000)
.panelPlot(seuratList, group_names, "Log10_nCount_RNA")
dev.off()
