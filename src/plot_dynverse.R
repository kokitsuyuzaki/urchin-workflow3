source("src/Functions3.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]
outfile4 <- commandArgs(trailingOnly=TRUE)[6]
# infile1 = 'output/hpbase/integrated/seurat_annotated.RData'
# infile2 = 'output/hpbase/integrated/dynverse/paga_celltype.RData'

# Load
load(infile1)
load(infile2)

# Plot
png(file = outfile1, width = 1200, height = 1200)
plot_dimred(model,
    grouping = seurat.integrated@meta.data$seurat_clusters)
dev.off()

png(file = outfile2, width = 1200, height = 1200)
plot_dimred(model,
    grouping = seurat.integrated@meta.data$sample)
dev.off()

png(file = outfile3, width = 1200, height = 1200)
plot_dimred(model,
    grouping = seurat.integrated@meta.data$celltype)
dev.off()

png(file = outfile4, width = 1200, height = 1200)
plot_dimred(model,
    grouping = seurat.integrated@meta.data$germlayer,
    color_density = "grouping")
dev.off()
