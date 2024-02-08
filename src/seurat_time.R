source("src/Functions.R")

# Parameter
set.seed(1234)
time <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

if(time == "24h"){
    infiles <- paste0("output/hpbase/",
        c("cont-24h", "DAPT-24h"), "/seurat.RData")
}
if(time == "36h"){
    infiles <- paste0("output/hpbase/",
        c("cont-36h", "DAPT-36h"), "/seurat.RData")
}
if(time == "48h"){
    infiles <- paste0("output/hpbase/",
        c("cont-48h", "DAPT-48h"), "/seurat.RData")
}
if(time == "72h"){
    infiles <- paste0("output/hpbase/",
        c("cont-72h", "DAPT-72h"), "/seurat.RData")
}
if(time == "96h"){
    infiles <- paste0("output/hpbase/",
        c("cont-96h", "DAPT-96h"), "/seurat.RData")
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
