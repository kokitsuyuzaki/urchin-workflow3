source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_integrated <- unlist(read.delim('plot/hpbase/integrated/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_integrated <- unlist(read.delim('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
load('output/hpbase/integrated/seurat.RData')

# Sort
names(all_states_integrated) <- NULL
names(bin_data_integrated) <- NULL

# 各データごとの状態No
target_integrated <- sapply(bin_data_integrated, function(x){
     which(all_states_integrated == x)
})

# Assign Labels
seurat.integrated$states <- target_integrated

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", label=FALSE, pt.size=1, label.size=6) + NoLegend()
ggsave(file=outfile1, g, dpi=200, width=6, height=6)

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", split.by="sample", label=FALSE, ncol=5, pt.size=1, label.size=6) + NoLegend()
ggsave(file=outfile2, g, dpi=200, width=20, height=12)
