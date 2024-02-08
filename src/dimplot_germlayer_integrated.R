source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# 色
cols <- sapply(names(table(seurat.integrated$germlayer)), function(x){
    target <- which(seurat.integrated$germlayer == x)[1]
    tmp <- seurat.integrated$germlayer_colors[target]
    names(tmp) <- seurat.integrated$germlayer[x]
    tmp
})
names(cols) <- gsub("\\.NA", "", names(cols))

# Idents切り替え
Idents(seurat.integrated) <- factor(seurat.integrated$germlayer, levels=names(cols))

# Plot
png(file=outfile1, width=1200, height=600)
p1 <- DimPlot(seurat.integrated, reduction = "umap", group.by="sample", label=TRUE, pt.size=2, label.size=6) + NoLegend()
p2 <- DimPlot(seurat.integrated, reduction = "umap", label=FALSE, pt.size=2, label.size=6, cols=cols)
p1 + p2
dev.off()

png(file=outfile2, width=2000, height=1000)
DimPlot(seurat.integrated, reduction = "umap", split.by="sample",
    ncol=5, label=FALSE, pt.size=2, label.size=6, cols=cols)
dev.off()
