source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile1)
load(infile2)

# Calculate the Percentage
seurat.integrated[["DoubletScore"]] <- dbl.dens

# Plot
png(file=outfile1, width=600, height=600)
FeaturePlot(seurat.integrated, features="DoubletScore",
    reduction = "umap", pt.size=2, label.size=6) + xlim(c(-15,15)) + ylim(c(-15,15))
dev.off()

seuratList <- .stratifySeurat(seurat.integrated, group_names)
png(file=outfile2, width=2000, height=1000)
.panelPlot(seuratList, group_names, "DoubletScore")
dev.off()

