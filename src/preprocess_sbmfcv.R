source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]

outfile1 <- commandArgs(trailingOnly=TRUE)[4]
outfile2 <- commandArgs(trailingOnly=TRUE)[5]
outfile3 <- commandArgs(trailingOnly=TRUE)[6]
outfile4 <- commandArgs(trailingOnly=TRUE)[7]
outfile5 <- commandArgs(trailingOnly=TRUE)[8]
outfile6 <- commandArgs(trailingOnly=TRUE)[9]
outfile7 <- commandArgs(trailingOnly=TRUE)[10]
outfile8 <- commandArgs(trailingOnly=TRUE)[11]
outfile9 <- commandArgs(trailingOnly=TRUE)[12]
outfile10 <- commandArgs(trailingOnly=TRUE)[13]
outfile11 <- commandArgs(trailingOnly=TRUE)[14]
outfile12 <- commandArgs(trailingOnly=TRUE)[15]
outfile13 <- commandArgs(trailingOnly=TRUE)[16]
outfile14 <- commandArgs(trailingOnly=TRUE)[17]
outfile15 <- commandArgs(trailingOnly=TRUE)[18]
outfile16 <- commandArgs(trailingOnly=TRUE)[19]
outfile17 <- commandArgs(trailingOnly=TRUE)[20]
outfile18 <- commandArgs(trailingOnly=TRUE)[21]
outfile19 <- commandArgs(trailingOnly=TRUE)[22]
outfile20 <- commandArgs(trailingOnly=TRUE)[23]
outfile21 <- commandArgs(trailingOnly=TRUE)[24]

# Loading
load(infile1)

# Preprocess
## UMAP
# ## PCA（だめだったやり方）
# out <- seurat.integrated@reductions$pca@cell.embeddings[,1:20]
out <- RunUMAP(seurat.integrated, dims=1:50, n.components=20L)@reductions$umap@cell.embeddings

## Normalization
# Min-Max Normalization（だめだったやり方）
# out <- out - min(out)
out <- apply(out, 2, function(x){x - min(x)})

# Cluster Number
# group <- as.character(Idents(seurat.integrated))
group <- as.character(seurat.integrated@meta.data$celltype)

# Time
time <- unlist(.timevec[seurat.integrated@meta.data$sample])

# Output
file.copy(infile1, outfile1)
file.copy(infile2, outfile2)
file.copy(infile3, outfile3)

write.table(out, outfile4, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out, outfile5, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out[grep("cont-", seurat.integrated@meta.data$sample), ],
    outfile6, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out[grep("cont-", seurat.integrated@meta.data$sample), ],
    outfile7, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out[grep("DAPT-", seurat.integrated@meta.data$sample), ],
    outfile8, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(out[grep("DAPT-", seurat.integrated@meta.data$sample), ],
    outfile9, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(group, outfile10, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group, outfile11, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile12, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile13, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile14, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(group[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile15, row.names=FALSE, col.names=FALSE, quote=FALSE)

write.table(time, outfile16, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time, outfile17, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile18, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("cont-", seurat.integrated@meta.data$sample)],
    outfile19, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile20, row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(time[grep("DAPT-", seurat.integrated@meta.data$sample)],
    outfile21, row.names=FALSE, col.names=FALSE, quote=FALSE)
