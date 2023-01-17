source("src/Functions.R")

# Parameter
set.seed(1234)
db <- commandArgs(trailingOnly=TRUE)[1]
sample <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]
indir <- paste0("output/", db, "/", sample, "/outs/filtered_feature_bc_matrix")

# Loading
data <- Read10X(data.dir=indir)

# Seurat Object
seurat.obj <- CreateSeuratObject(counts = data)

# Normalization
seurat.obj <- NormalizeData(seurat.obj)
seurat.obj <- FindVariableFeatures(seurat.obj)
seurat.obj <- ScaleData(seurat.obj)

# Dimensional Reduction
seurat.obj <- RunPCA(seurat.obj)
seurat.obj <- RunUMAP(seurat.obj, dims=1:30)

# Clustering
seurat.obj <- FindNeighbors(seurat.obj, dims=1:30)
seurat.obj <- FindClusters(seurat.obj)

# Save
save(seurat.obj, file=outfile)