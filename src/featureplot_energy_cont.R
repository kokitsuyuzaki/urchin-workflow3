source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_cont <- unlist(read.delim('plot/hpbase/cont/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/hpbase/cont/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_cont <- unlist(read.table('plot/hpbase/cont/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/cont/seurat.RData')

# Sort
names(all_states_cont) <- NULL
names(bin_data_cont) <- NULL

target_cont <- sapply(bin_data_cont, function(x){
     which(all_states_cont == x)
})

energy_cont_sorted <- energy_cont[target_cont]

# Assign Labels
seurat.integrated$energy <- energy_cont_sorted

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=2)
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", split.by="sample", ncol=5, pt.size=2)
png(file=outfile2, width=3000, height=600)
print(g)
dev.off()

