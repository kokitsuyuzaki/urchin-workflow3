source("src/Functions2.R")

# Parameter
db <- commandArgs(trailingOnly=TRUE)[1]
infile <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)
infile2 <- paste0("output/", db, "/integrated/seurat.RData")
load(infile2)

# Embed Velocity on UMAP coordinates
umap.coordinates <- seurat.integrated@reductions$umap@cell.embeddings
rownames(umap.coordinates) <- colnames(velo.obj)
colnames(umap.coordinates) <- c("V1", "V2")
reducedDims(velo.obj) <- SimpleList(UMAP = umap.coordinates)

# Embedding
em <- embedVelocity(reducedDim(velo.obj, "UMAP"), velo.obj)

if(length(which(is.nan(em))) == 0){
    # Plot
    png(file=outfile, width=600, height=600)
    plot(plotVelocityStream(velo.obj, em, color.streamlines=TRUE, use.dimred="UMAP"))
    dev.off()
}else{
    file.create(outfile)
}