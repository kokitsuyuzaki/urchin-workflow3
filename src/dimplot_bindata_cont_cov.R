source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
bindata_cont <- read.table('output/hpbase/cont/sbmfcv/BIN_DATA.tsv', header=FALSE)
load('output/hpbase/cont_stratified/seurat.RData')

# Assign Labels
for(i in seq_len(ncol(bindata_cont))){
     cmd <- paste0("seurat.integrated$bindata_", i, " <- bindata_cont[, i]")
     eval(parse(text=cmd))
}

# Setting
dir.create("plot/hpbase/cont_cov/bindata/")

# Plot cont
for(i in seq_len(ncol(bindata_cont))){
     filename1 <- paste0("plot/hpbase/cont_cov/bindata/", i, ".png")
     filename2 <- paste0("plot/hpbase/cont_cov/bindata/", i, "_splitby.png")
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
