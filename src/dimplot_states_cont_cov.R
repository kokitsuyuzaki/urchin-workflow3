source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_cont <- unlist(read.delim('plot/hpbase/cont_cov/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/hpbase/cont/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
load('output/hpbase/cont_stratified/seurat_annotated.RData')

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
seurat.integrated <- seurat.integrated[, intersect(target1, target2)]

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
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile1, width=600, height=600)
print(g)
dev.off()

# Plot
g <- DimPlot(seurat.integrated, reduction = "umap", group.by="states", split.by="sample", label=FALSE, pt.size=2, label.size=6) + NoLegend()
png(file=outfile2, width=2400, height=600)
print(g)
dev.off()
