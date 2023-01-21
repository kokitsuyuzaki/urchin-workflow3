source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Marker detection
seurat_clusters <- seurat.integrated@meta.data$seurat_clusters
outList <- lapply(seq(nlevels(seurat_clusters)), function(i){
     tmp <- FindConservedMarkers(seurat.integrated, ident.1=levels(seurat_clusters)[i],
         grouping.var = "sample")
     genenames <- rownames(tmp)
     cbind(genenames, tmp)
})
names(outList) <- paste0("cluster", seq(nlevels(seurat_clusters)))

# Save
write_xlsx(outList, outfile)
