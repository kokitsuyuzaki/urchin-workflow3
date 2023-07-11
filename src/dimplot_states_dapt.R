source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_dapt <- unlist(read.delim('plot/hpbase/dapt/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_dapt <- unlist(read.delim('output/hpbase/dapt/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
load('output/hpbase/dapt/seurat.RData')

# Sort
names(all_states_dapt) <- NULL
names(bin_data_dapt) <- NULL

# 各データごとの状態No
target_dapt <- sapply(bin_data_dapt, function(x){
     which(all_states_dapt == x)
})

# Assign Labels
seurat.dapt$states <- target_dapt

# Plot
g <- DimPlot(seurat.dapt, reduction = "umap", group.by="states", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.dapt, reduction = "umap", group.by="states", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile2, width=2400, height=600)
print(g)
dev.off()

