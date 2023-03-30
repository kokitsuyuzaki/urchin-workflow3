source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# Plot
png(file=outfile1, width=1200, height=600)
p1 <- DimPlot(seurat.integrated, reduction = "umap", group.by="sample", label=TRUE, pt.size=2, label.size=6) + NoLegend()
p2 <- DimPlot(seurat.integrated, reduction = "umap", label=FALSE, pt.size=2, label.size=6, cols=cols)
p1 + p2
dev.off()

png(file=outfile2, width=2000, height=1000)
DimPlot(seurat.integrated, reduction = "umap", split.by="sample",
    ncol=5, label=FALSE, pt.size=2, label.size=6, cols=cols)
dev.off()
