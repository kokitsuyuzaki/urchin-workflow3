source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Plot
png(file=outfile, width=600, height=600)
DotPlot(seurat.obj, features = markers) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
