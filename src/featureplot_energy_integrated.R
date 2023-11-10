source("src/Functions.R")

# Parameter
outfile1 <- commandArgs(trailingOnly=TRUE)[1]
outfile2 <- commandArgs(trailingOnly=TRUE)[2]

# Loading
all_states_integrated <- unlist(read.delim('plot/hpbase/integrated/Landscaper/Allstates.tsv', header=FALSE, sep="|"))
bin_data_integrated <- unlist(read.delim('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE, sep="|"))
energy_integrated <- unlist(read.table('plot/hpbase/integrated/Landscaper/E.tsv', header=FALSE))
load('output/hpbase/integrated/seurat.RData')

# rm \t
bin_data_integrated <- gsub("\t", " ", bin_data_integrated)

# Sort
names(all_states_integrated) <- NULL
names(bin_data_integrated) <- NULL

target_integrated <- sapply(bin_data_integrated, function(x){
     which(all_states_integrated == x)
})

energy_integrated_sorted <- energy_integrated[target_integrated]

# Assign Labels
seurat.integrated$energy <- energy_integrated_sorted

# Plot
g <- FeaturePlot(seurat.integrated, reduction = "umap", features="energy", pt.size=1)
ggsave(file=outfile1, g, dpi=200, width=6, height=6)

# Plot
seuratList <- .stratifySeurat(seurat.integrated, group_names)
g <- .panelPlot(seuratList, group_names, "energy")
ggsave(file=outfile2, g, dpi=200, width=30, height=15)
