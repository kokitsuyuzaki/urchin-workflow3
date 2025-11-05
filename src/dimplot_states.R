source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile4 <- commandArgs(trailingOnly=TRUE)[4]

# Loading
all_states_cont <- unlist(read.delim('plot/hpbase/cont/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
all_states_dapt <- unlist(read.delim('plot/hpbase/dapt/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/hpbase/cont/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
bin_data_dapt <- unlist(read.delim('output/hpbase/dapt/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
load('output/hpbase/cont/seurat.RData')
load('output/hpbase/dapt/seurat.RData')

# Sort
names(all_states_cont) <- NULL
names(all_states_dapt) <- NULL
names(bin_data_cont) <- NULL
names(bin_data_dapt) <- NULL

# 各データごとの状態No
target_cont <- sapply(bin_data_cont, function(x){
     which(all_states_cont == x)
})
target_dapt <- sapply(bin_data_dapt, function(x){
     which(all_states_dapt == x)
})

# Assign Labels
names(target_cont) <- colnames(seurat.cont)
names(target_dapt) <- colnames(seurat.dapt)
seurat.cont$states <- target_cont
seurat.dapt$states <- target_dapt

# Setting
dir.create("plot/hpbase/integrated/states/")

# Plot
g <- DimPlot(seurat.cont, reduction = "umap", group.by="states", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.dapt, reduction = "umap", group.by="states", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile2, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.cont, reduction = "umap", group.by="states", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile3, width=2000, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.dapt, reduction = "umap", group.by="states", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile4, width=2000, height=600)
print(g)
dev.off()
