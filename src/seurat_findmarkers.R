source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Marker detection
seurat_clusters <- seurat.obj@meta.data$seurat_clusters
outList <- lapply(seq(nlevels(seurat_clusters)), function(i){
     FindMarkers(seurat.obj, ident.1=levels(seurat_clusters)[i])
})
names(outList) <- paste0("cluster", seq(nlevels(seurat_clusters)))

# Save
write_xlsx(outList, outfile)
