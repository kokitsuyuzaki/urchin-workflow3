source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
bindata_integrated <- read.table('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE)
load('output/hpbase/integrated/seurat.RData')

# Assign Labels
for(i in seq_len(ncol(bindata_integrated))){
     cmd <- paste0("seurat.integrated$bindata_", i, " <- bindata_integrated[, i]")
     eval(parse(text=cmd))
}

# Setting
dir.create("plot/hpbase/integrated/bindata/")

# Plot integrated
for(i in seq_len(ncol(bindata_integrated))){
     filename1 <- paste0("plot/hpbase/integrated/bindata/", i, ".png")
     filename2 <- paste0("plot/hpbase/integrated/bindata/", i, "_splitby.png")
     groupname <- paste0("bindata_", i)
     # Plot
     g <- DimPlot(seurat.integrated, reduction = "umap", group.by=groupname, label=TRUE, pt.size=1, label.size=6, cols=c(4,2)) + NoLegend()
     ggsave(file=filename1, g, dpi=200, width=6, height=6)
     g <- DimPlot(seurat.integrated, reduction = "umap", group.by=groupname, split.by="sample",
         ncol=5, label=TRUE, pt.size=1, label.size=6, cols=c(4,2)) + NoLegend()
     ggsave(file=filename2, g, dpi=200, width=20, height=12)
}

# Save
file.create(outfile)
