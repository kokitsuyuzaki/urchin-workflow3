source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# 色
cols <- sapply(names(table(seurat.obj$germlayer)), function(x){
    target <- which(seurat.obj$germlayer == x)[1]
    tmp <- seurat.obj$germlayer_colors[target]
    names(tmp) <- seurat.obj$germlayer[x]
    tmp
})
names(cols) <- gsub("\\.NA", "", names(cols))

# Idents切り替え
Idents(seurat.obj) <- factor(seurat.obj$germlayer, levels=names(cols))

# Plot
png(file=outfile, width=600, height=600)
DimPlot(seurat.obj, label=FALSE, pt.size=2, label.size=6, cols=cols)
dev.off()
