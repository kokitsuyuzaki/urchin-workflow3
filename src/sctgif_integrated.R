source("src/Functions.R")

# Parameter
db <- commandArgs(trailingOnly=TRUE)[1]
go <- commandArgs(trailingOnly=TRUE)[2]
infile1 <- commandArgs(trailingOnly=TRUE)[3]
infile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile <- commandArgs(trailingOnly=TRUE)[5]
out.dir <- gsub("index.html", "", outfile)

# Loading
load(infile1)
load(infile2)
gmt <- eval(parse(text=paste0("gmt_go_", go, "_", db)))

# Setting
sce <- as.SingleCellExperiment(seurat.integrated)
settingTGIF(sce, gmt, reducedDimNames="UMAP", assayNames="logcounts", nbins=100)

# Calculation
calcTGIF(sce, ndim=5)

# Save
reportTGIF(sce, out.dir=out.dir,
    html.open=FALSE, title="scTGIF Report for sea urchin",
    author="Koki Tsuyuzaki")
