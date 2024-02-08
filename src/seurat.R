source("src/Functions.R")

# Parameter
set.seed(1234)
sample <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
indir <- paste0("output/hpbase/", sample, "/outs/filtered_feature_bc_matrix")

# Loading
geneid_to_genename <- read.csv('data/geneid_to_genename.csv', row.names=1)
data <- Read10X(data.dir=indir)

# Gene ID => Gene Name
tmp <- convertRowID(
    input = as.matrix(data),
    rowID = rownames(data),
    LtoR = geneid_to_genename[, 1:2])
data <- as(tmp$output, "sparseMatrix")

# ここで、"_"を"-"に勝手に変換するSeuratの仕様と組み合わさって、Featureが重複してコケる
data <- data[grep(blacklist_gene, rownames(data), invert=TRUE), ]

# Seurat Object
seurat.obj <- CreateSeuratObject(counts = data)

# SCTransform
seurat.obj <- SCTransform(seurat.obj)

# Dimensional Reduction
seurat.obj <- RunPCA(seurat.obj)
seurat.obj <- RunUMAP(seurat.obj, dims=1:30)

# Clustering
seurat.obj <- FindNeighbors(seurat.obj, dims=1:30)
seurat.obj <- FindClusters(seurat.obj)

# Save
save(seurat.obj, file=outfile)