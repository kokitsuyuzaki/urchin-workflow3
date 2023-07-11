source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_dapt <- unlist(read.delim('plot/hpbase/dapt/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_dapt <- unlist(read.delim('output/hpbase/dapt/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_dapt <- unlist(read.table('plot/hpbase/dapt/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/dapt/seurat.RData')

# Sort
names(all_states_dapt) <- NULL
names(bin_data_dapt) <- NULL

target_dapt <- sapply(bin_data_dapt, function(x){
     which(all_states_dapt == x)
})

energy_dapt_sorted <- energy_dapt[target_dapt]

# Assign Labels
seurat.dapt$energy <- energy_dapt_sorted

# Plot
g <- FeaturePlot(seurat.dapt, reduction = "umap", features="energy", pt.size=2)
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- FeaturePlot(seurat.dapt, reduction = "umap", features="energy", split.by="sample", ncol=5, pt.size=2)
png(file=outfile2, width=3000, height=600)
print(g)
dev.off()

