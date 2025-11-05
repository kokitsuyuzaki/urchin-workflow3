source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
bindata_DAPT <- read.table('output/hpbase/DAPT/sbmfcv/BIN_DATA.tsv', header=FALSE)
load('output/hpbase/DAPT_stratified/seurat_annotated.RData')

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
seurat.integrated <- seurat.integrated[, intersect(target1, target2)]

# Assign Labels
for(i in seq_len(ncol(bindata_DAPT))){
     cmd <- paste0("seurat.integrated$bindata_", i, " <- bindata_DAPT[, i]")
     eval(parse(text=cmd))
}

# Setting
dir.create("plot/hpbase/DAPT/bindata/")

# Plot DAPT
for(i in seq_len(ncol(bindata_DAPT))){
     filename1 <- paste0("plot/hpbase/DAPT/bindata/", i, ".png")
     filename2 <- paste0("plot/hpbase/DAPT/bindata/", i, "_splitby.png")
     groupname <- paste0("bindata_", i)
     # Plot
     g <- DimPlot(seurat.integrated, reduction = "umap", group.by=groupname, label=TRUE, pt.size=2, label.size=6, cols=c(4,2)) + NoLegend()
     png(file=filename1, width=600, height=600)
     print(g)
     dev.off()
     # Plot
     g <- DimPlot(seurat.integrated, reduction = "umap", group.by=groupname, split.by="sample",
         ncol=5, label=TRUE, pt.size=2, label.size=6, cols=c(4,2)) + NoLegend()
     png(file=filename2, width=2400, height=600)
     print(g)
     dev.off()
}

# Save
file.create(outfile)
