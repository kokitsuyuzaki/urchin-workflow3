source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
bindata_dapt <- read.table('output/hpbase/dapt/sbmfcv/BIN_DATA.tsv', header=FALSE)
load('output/hpbase/dapt/seurat.RData')

# Assign Labels
for(i in seq_len(ncol(bindata_dapt))){
     cmd <- paste0("seurat.dapt$bindata_", i, " <- bindata_dapt[, i]")
     eval(parse(text=cmd))
}

# Setting
dir.create("plot/hpbase/dapt/bindata/")

# Plot dapt
for(i in seq_len(ncol(bindata_dapt))){
     filename1 <- paste0("plot/hpbase/dapt/bindata/", i, ".png")
     filename2 <- paste0("plot/hpbase/dapt/bindata/", i, "_splitby.png")
     groupname <- paste0("bindata_", i)
     # Plot
     g <- DimPlot(seurat.dapt, reduction = "umap", group.by=groupname, label=TRUE, pt.size=2, label.size=6, cols=c(4,2)) + NoLegend()
     png(file=filename1, width=600, height=600)
     print(g)
     dev.off()
     # Plot
     g <- DimPlot(seurat.dapt, reduction = "umap", group.by=groupname, split.by="sample",
         ncol=5, label=TRUE, pt.size=2, label.size=6, cols=c(4,2)) + NoLegend()
     png(file=filename2, width=2400, height=600)
     print(g)
     dev.off()
}

# Save
file.create(outfile)
