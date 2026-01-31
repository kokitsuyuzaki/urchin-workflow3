source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
all_states_cont <- unlist(read.delim('plot/echinobase/Oulhen/cont/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/echinobase/Oulhen/cont/binpca/BIN_DATA.tsv', header=FALSE, sep="|"))
load('output/echinobase/Oulhen/cont/seurat_annotated_landscaper.RData')

# Sort
names(all_states_cont) <- NULL
names(bin_data_cont) <- NULL

# 各データごとの状態No
target_cont <- sapply(bin_data_cont, function(x){
     which(all_states_cont == x)
})

# Assign Labels
names(target_cont) <- colnames(seurat.integrated)
seurat.integrated$states <- target_cont

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", label=FALSE, pt.size=3, label.size=6) + NoLegend() +
theme(axis.line = element_blank(),
           axis.text.x = element_blank(),
           axis.text.y = element_blank(),
           axis.ticks = element_blank(),
           axis.title.x = element_blank(),
           axis.title.y = element_blank(),
           panel.background = element_blank(),
           panel.border = element_blank(),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.background = element_blank())
png(file=outfile, width=600, height=600)
print(g)
dev.off()
