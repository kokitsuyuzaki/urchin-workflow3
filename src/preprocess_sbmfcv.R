source("src/Functions.R")

# Parameter
infile <- 'output/hpbase/integrated/seurat_annotated.RData'
outfile1 <- 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
outfile2 <- 'output/hpbase/integrated/seurat.tsv'
outfile3 <- 'output/hpbase/integrated/group.tsv'
outfile4 <- 'output/hpbase/integrated_cov/group.tsv'
outfile5 <- 'output/hpbase/cont/group.tsv'
outfile6 <- 'output/hpbase/cont_cov/group.tsv'
outfile7 <- 'output/hpbase/DAPT/group.tsv'
outfile8 <- 'output/hpbase/DAPT_cov/group.tsv'
outfile9 <- 'output/hpbase/integrated/cov.tsv'
outfile10 <- 'output/hpbase/integrated_cov/cov.tsv'
outfile11 <- 'output/hpbase/cont/cov.tsv'
outfile12 <- 'output/hpbase/cont_cov/cov.tsv'
outfile13 <- 'output/hpbase/DAPT/cov.tsv'
outfile14 <- 'output/hpbase/DAPT_cov/cov.tsv'

# Loading
load(infile)

# Preprocess
## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
target <- intersect(target1, target2)
seurat.integrated <- seurat.integrated[, target]

# Dimensional Reduction
## 純粋にNMFのみ（だめだったやり方）
# out <- as.matrix(seurat.integrated@assays$RNA@data[VariableFeatures(seurat.integrated), ])

# ## PCA（だめだったやり方）
# out <- seurat.integrated@reductions$pca@cell.embeddings[,1:20]

## PCA & UMAPやり直し
seurat.integrated <- RunPCA(seurat.integrated)
seurat.integrated <- RunUMAP(seurat.integrated, dims=1:30)
out <- seurat.integrated@reductions$pca@cell.embeddings[, 1:10]

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
save(seurat.integrated, file=outfile1)
write.table(out, outfile2, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group, outfile3, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group, outfile4, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile5, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile6, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile7, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile8, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(time, outfile9, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time, outfile10, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile11, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile12, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile13, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile14, row.names=FALSE, col.names=FALSE, quote=FALSE)
