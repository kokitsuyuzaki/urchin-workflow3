source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_dapt <- unlist(read.delim('plot/hpbase/dapt/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_dapt <- unlist(read.delim('output/hpbase/dapt/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
basin_dapt <- unlist(read.table('plot/hpbase/dapt/Landscaper/Basin.tsv', header=FALSE))
load('output/hpbase/dapt/seurat.RData')

# Sort
names(all_states_dapt) <- NULL
names(bin_data_dapt) <- NULL

# 各データごとの状態No
target_dapt <- sapply(bin_data_dapt, function(x){
     which(all_states_dapt == x)
})

basin_dapt_sorted <- rep(0, length=length(bin_data_dapt))

target_dapt2 <- unlist(sapply(basin_dapt, function(x){
     which(target_dapt == x)
}))

basin_dapt_sorted[target_dapt2] <- 1

# Assign Labels
seurat.dapt$basin <- basin_dapt_sorted

# Plot
g <- DimPlot(seurat.dapt, reduction = "umap", group.by="basin", label=FALSE, pt.size=2, label.size=6, cols=c(8,3)) + NoLegend()
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.dapt, reduction = "umap", group.by="basin", split.by="sample", label=FALSE, pt.size=2, label.size=6, cols=c(8,3)) + NoLegend()
png(file=outfile2, width=2400, height=600)
print(g)
dev.off()
