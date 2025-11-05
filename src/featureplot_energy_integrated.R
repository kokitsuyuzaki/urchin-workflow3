source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile4 <- commandArgs(trailingOnly=TRUE)[4]
outfile5 <- commandArgs(trailingOnly=TRUE)[5]
outfile6 <- commandArgs(trailingOnly=TRUE)[6]

####################### DAPT Step ####################
all_states_DAPT <- unlist(read.delim('plot/hpbase/integrated/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_DAPT <- unlist(read.delim('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_DAPT <- unlist(read.table('plot/hpbase/integrated/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/DAPT_stratified/seurat_annotated.RData')

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
target <- intersect(target1, target2)
seurat.integrated <- seurat.integrated[, target]

bin_data_DAPT <- gsub("\t", " ", bin_data_DAPT)
bin_data_DAPT <- bin_data_DAPT[target]

# Sort
names(all_states_DAPT) <- NULL
names(bin_data_DAPT) <- NULL

target_DAPT <- sapply(bin_data_DAPT, function(x){
     which(all_states_DAPT == x)
})

energy_DAPT_sorted <- energy_DAPT[target_DAPT]
energy_DAPT_sorted <- as.matrix(energy_DAPT_sorted)
rownames(energy_DAPT_sorted) <- colnames(seurat.integrated)
######################################################

####################### cont Step ####################
all_states_cont <- unlist(read.delim('plot/hpbase/integrated/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_cont <- unlist(read.delim('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_cont <- unlist(read.table('plot/hpbase/integrated/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/cont_stratified/seurat_annotated.RData')

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
target <- intersect(target1, target2)
seurat.integrated <- seurat.integrated[, target]

bin_data_cont <- gsub("\t", " ", bin_data_cont)
bin_data_cont <- bin_data_cont[target]

# Sort
names(all_states_cont) <- NULL
names(bin_data_cont) <- NULL

target_cont <- sapply(bin_data_cont, function(x){
     which(all_states_cont == x)
})

energy_cont_sorted <- energy_cont[target_cont]
energy_cont_sorted <- as.matrix(energy_cont_sorted)
rownames(energy_cont_sorted) <- colnames(seurat.integrated)
######################################################

####################### integrated Step ####################
all_states_integrated <- unlist(read.delim('plot/hpbase/integrated/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_integrated <- unlist(read.delim('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_integrated <- unlist(read.table('plot/hpbase/integrated/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/integrated/seurat_annotated.RData')

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
target <- intersect(target1, target2)
seurat.integrated <- seurat.integrated[, target]

bin_data_integrated <- gsub("\t", " ", bin_data_integrated)
bin_data_integrated <- bin_data_integrated[target]

# Sort
names(all_states_integrated) <- NULL
names(bin_data_integrated) <- NULL

target_integrated <- sapply(bin_data_integrated, function(x){
     which(all_states_integrated == x)[1]
})

energy_integrated_sorted <- energy_integrated[target_integrated]
energy_integrated_sorted <- as.matrix(energy_integrated_sorted)
rownames(energy_integrated_sorted) <- colnames(seurat.integrated)
######################################################

# Assign Labels
seurat.integrated$energy <- energy_integrated_sorted
seurat.integrated$cont_DAPT_energy <- 0
seurat.integrated$cont_DAPT_energy[rownames(energy_cont_sorted)] <- energy_cont_sorted
seurat.integrated$cont_DAPT_energy[rownames(energy_DAPT_sorted)] <- energy_DAPT_sorted

# Seurat => SingleCellExperiment
sce <- SingleCellExperiment(assays = list(counts = seurat.integrated@assays$RNA@counts))
reducedDims(sce) <- list(UMAP = seurat.integrated@reductions$umap@cell.embeddings)
logcounts(sce) <- log10(counts(sce) + 1)

# Assign Condition Labels
colData(sce)$condition <- gsub("-.*", "", seurat.integrated@meta.data$sample)

# Assign Energy
colData(sce)$energy <- seurat.integrated$energy
colData(sce)$cont_DAPT_energy <- seurat.integrated$cont_DAPT_energy

# Setting Hex bin
sce <- make_hexbin(sce, nbins = 30, dimension_reduction = "UMAP")
save(sce, file=outfile1)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100))
ggsave(file=outfile2, g, dpi=200, width=6, height=6)

# Plot
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "")
ggsave(file=outfile3, g, dpi=200, width=6, height=6)

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1) + scale_colour_gradientn(colours = viridis(100), limits=c(-7.0, 8.5))
ggsave(file=outfile4, g, dpi=200, width=6, height=6)

# Plot
g <- plot_hexbin_meta(sce, col = "energy", action = "median", title = "")
ggsave(file=outfile5, g, dpi=200, width=6, height=6)

# Plot
seuratList <- .stratifySeurat(seurat.integrated, group_names)
g <- .panelPlotMeta(seuratList, group_names, "energy")
ggsave(file=outfile6, g, dpi=200, width=30, height=15)
