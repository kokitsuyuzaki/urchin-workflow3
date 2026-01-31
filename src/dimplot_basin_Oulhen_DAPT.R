source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
all_states_DAPT <- unlist(read.delim('plot/echinobase/Oulhen/DAPT/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_DAPT <- unlist(read.delim('output/echinobase/Oulhen/DAPT/binpca/BIN_DATA.tsv', header=FALSE, sep="|"))
basin_DAPT <- unlist(read.table('plot/echinobase/Oulhen/DAPT/Landscaper/Basin.tsv', header=FALSE))
load('output/echinobase/Oulhen/DAPT/seurat_annotated_landscaper.RData')

# Sort
names(all_states_DAPT) <- NULL
names(bin_data_DAPT) <- NULL

# 各データごとの状態No
target_DAPT <- sapply(bin_data_DAPT, function(x){
     which(all_states_DAPT == x)
})

basin_DAPT_sorted <- rep(0, length=length(bin_data_DAPT))

target_DAPT2 <- unlist(sapply(basin_DAPT, function(x){
     which(target_DAPT == x)
}))

basin_DAPT_sorted[target_DAPT2] <- 1

# Assign Labels
seurat.integrated$basin <- basin_DAPT_sorted

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="basin", label=FALSE, pt.size=3, label.size=6, cols=c(8,3)) + NoLegend() +
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
