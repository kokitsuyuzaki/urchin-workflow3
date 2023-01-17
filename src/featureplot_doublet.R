source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
load(infile2)

# Plot
seurat.obj[["DoubletScore"]] <- dbl.dens

# Plot
png(file=outfile, width=600, height=600)
FeaturePlot(seurat.obj, features="DoubletScore", pt.size=2, label.size=6)
dev.off()
