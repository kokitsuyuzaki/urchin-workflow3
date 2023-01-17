source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# Plot
png(file=outfile1, width=600, height=600)
DotPlot(seurat.integrated, features = markers) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

png(file=outfile2, width=1000, height=2000)
DotPlot(seurat.integrated, features = markers, split.by="sample", cols=sample_colors) + theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()
