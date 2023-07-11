source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_cont <- unlist(read.delim('plot/hpbase/cont/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/hpbase/cont/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
load('output/hpbase/cont/seurat.RData')

# Sort
names(all_states_cont) <- NULL
names(bin_data_cont) <- NULL

# 各データごとの状態No
target_cont <- sapply(bin_data_cont, function(x){
     which(all_states_cont == x)
})

# Assign Labels
seurat.cont$states <- target_cont

# Plot
g <- DimPlot(seurat.cont, reduction = "umap", group.by="states", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.cont, reduction = "umap", group.by="states", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile2, width=2400, height=600)
print(g)
dev.off()
