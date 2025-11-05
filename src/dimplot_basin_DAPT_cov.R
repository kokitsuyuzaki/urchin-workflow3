source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_DAPT <- unlist(read.delim('plot/hpbase/DAPT_cov/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_DAPT <- unlist(read.delim('output/hpbase/DAPT/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
basin_DAPT <- unlist(read.table('plot/hpbase/DAPT_cov/Landscaper/Basin.tsv', header=FALSE))
load('output/hpbase/DAPT_stratified/seurat_annotated.RData')

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
seurat.integrated <- seurat.integrated[, intersect(target1, target2)]

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
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="basin", label=FALSE, pt.size=2, label.size=6, cols=c(8,3)) + NoLegend()
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="basin", split.by="sample", label=FALSE, pt.size=2, label.size=6, cols=c(8,3)) + NoLegend()
png(file=outfile2, width=2400, height=600)
print(g)
dev.off()
