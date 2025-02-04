source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile3 <- commandArgs(trailingOnly=TRUE)[3]

# Loading
all_states_DAPT <- unlist(read.delim('plot/hpbase/DAPT_cov/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_DAPT <- unlist(read.delim('output/hpbase/DAPT/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_DAPT <- unlist(read.table('plot/hpbase/DAPT_cov/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/DAPT_stratified/seurat.RData')

# Sort
names(all_states_DAPT) <- NULL
names(bin_data_DAPT) <- NULL

target_DAPT <- sapply(bin_data_DAPT, function(x){
     which(all_states_DAPT == x)
})

energy_DAPT_sorted <- energy_DAPT[target_DAPT]

# Assign Labels
seurat.integrated$energy <- energy_DAPT_sorted

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=2) + scale_colour_gradientn(colours = viridis(100))
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=2) + scale_colour_gradientn(colours = viridis(100), limits=c(-7.0, 8.5))
png(file=outfile2, width=600, height=600)
print(g)
dev.off()

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", split.by="sample", ncol=5, pt.size=2)
png(file=outfile3, width=3000, height=600)
print(g)
dev.off()

