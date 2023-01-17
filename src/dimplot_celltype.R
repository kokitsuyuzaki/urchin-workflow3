source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Plot
png(file=outfile, width=600, height=600)
DimPlot(seurat.obj, label=FALSE, pt.size=2, label.size=6, cols=cols)
dev.off()
