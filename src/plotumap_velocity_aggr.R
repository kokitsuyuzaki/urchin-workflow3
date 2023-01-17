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
grid.df <- gridVectors(reducedDim(velo.obj, "UMAP"), em)

if(length(which(is.nan(em))) == 0){
	# Plot
	png(file=outfile, width=600, height=600)
	plot(plotUMAP(velo.obj, colour_by="velocity_pseudotime") +
	    geom_segment(data=grid.df, mapping=aes(x=start.V1, y=start.V2,
	        xend=end.V1, yend=end.V2), arrow=arrow(length=unit(0.1, "inches"))))
	dev.off()
}else{
	file.create(outfile)
}
