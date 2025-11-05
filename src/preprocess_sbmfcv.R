source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]
outfile4 <- commandArgs(trailingOnly=TRUE)[5]
outfile5 <- commandArgs(trailingOnly=TRUE)[6]
outfile6 <- commandArgs(trailingOnly=TRUE)[7]
outfile7 <- commandArgs(trailingOnly=TRUE)[8]
outfile8 <- commandArgs(trailingOnly=TRUE)[9]
outfile9 <- commandArgs(trailingOnly=TRUE)[10]
outfile10 <- commandArgs(trailingOnly=TRUE)[11]
outfile11 <- commandArgs(trailingOnly=TRUE)[12]
outfile12 <- commandArgs(trailingOnly=TRUE)[13]
outfile13 <- commandArgs(trailingOnly=TRUE)[14]

# Loading
load(infile1)

# Preprocess
## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
seurat.integrated <- seurat.integrated[, intersect(target1, target2)]

# Dimensional Reduction
## 純粋にNMFのみ（だめだったやり方）
# out <- as.matrix(seurat.integrated@assays$RNA@data[VariableFeatures(seurat.integrated), ])
# print(dim(out))

# ## PCA（だめだったやり方）
# out <- seurat.integrated@reductions$pca@cell.embeddings[,1:20]

## UMAP
## PCAやり直し
seurat.integrated <- RunPCA(seurat.integrated)
out <- seurat.integrated@reductions$pca@cell.embeddings[,1:10]

# set.seed(123456)
# out <- RunUMAP(seurat.integrated, dims=1:50, n.components=20L)@reductions$umap@cell.embeddings

# ## Normalization
# # Min-Max Normalization（だめだったやり方）
# # out <- out - min(out)
# out <- apply(out, 2, function(x){x - min(x)})

# Cluster Number
# group <- as.character(Idents(seurat.integrated))
group <- as.character(seurat.integrated@meta.data$celltype)

# Time
time <- unlist(.timevec[seurat.integrated@meta.data$sample])

# Output
write.table(out, outfile1, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group, outfile2, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group, outfile3, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile4, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile5, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile6, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile7, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(time, outfile8, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time, outfile9, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile10, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile11, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile12, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile13, row.names=FALSE, col.names=FALSE, quote=FALSE)
