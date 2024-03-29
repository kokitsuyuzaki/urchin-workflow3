source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# 色
cols <- sapply(names(table(seurat.obj$celltype)), function(x){
    target <- which(seurat.obj$celltype == x)[1]
    tmp <- seurat.obj$celltype_colors[target]
    names(tmp) <- seurat.obj$celltype[x]
    tmp
})
names(cols) <- gsub("\\.NA", "", names(cols))

# Idents切り替え
Idents(seurat.obj) <- factor(seurat.obj$celltype, levels=names(cols))

# Plot
png(file=outfile, width=600, height=600)
DimPlot(seurat.obj, label=FALSE, pt.size=2, label.size=6, cols=cols)
dev.off()
