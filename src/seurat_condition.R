source("src/Functions.R")

# Parameter
set.seed(1234)
condition <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

if(condition == "cont"){
    infiles <- paste0("output/hpbase/",
        c("cont-24h", "cont-36h", "cont-48h", "cont-72h", "cont-96h"),
        "/seurat.RData")
}else{
    infiles <- paste0("output/hpbase/",
        c("DAPT-24h", "DAPT-36h", "DAPT-48h", "DAPT-72h", "DAPT-96h"),
        "/seurat.RData")
}

# Loading
seurat.list <- lapply(infiles, function(x){
    load(x)
    return(seurat.obj)
})

# Integration
features <- SelectIntegrationFeatures(object.list = seurat.list)
anchors <- FindIntegrationAnchors(object.list = seurat.list, anchor.features = features)
seurat.integrated <- IntegrateData(anchorset = anchors)
DefaultAssay(seurat.integrated) <- "integrated"

# SCTransform
seurat.integrated <- SCTransform(seurat.integrated)

# Dimension Reduction
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, dims=1:30)

# Clustering
seurat.integrated <- FindNeighbors(seurat.integrated, dims=1:30)
seurat.integrated <- FindClusters(seurat.integrated)

# Group
seurat.integrated@meta.data$sample <- unlist(lapply(seq_along(seurat.list),
    function(x){
        rep(group_names[x], length=ncol(seurat.list[[x]]))
    }))

# Save
save(seurat.integrated, file=outfile)
