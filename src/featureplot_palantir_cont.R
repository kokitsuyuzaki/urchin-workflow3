source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile1)
pseudotime <- read.csv(infile2, header=TRUE, row.names=1)

# Assign Labels
seurat.integrated[["pseudotime"]] <- pseudotime

# Plot
png(file=outfile1, width=600, height=600)
FeaturePlot(seurat.integrated, features="pseudotime",
    reduction = "umap", pt.size=2, label.size=6) + xlim(c(-15,15)) + ylim(c(-15,15)) + scale_color_viridis(option = "inferno", limits = c(0,1)) 
dev.off()

seuratList <- .stratifySeurat2(seurat.integrated, c("cont-36h", "cont-48h", "cont-72h", "cont-96h"))
png(file=outfile2, width=2000, height=650)
.panelPlot2(seuratList, c("cont-36h", "cont-48h", "cont-72h", "cont-96h"), "pseudotime")
dev.off()
