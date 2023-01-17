source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Predict Doublets
cds <- as.cell_data_set(seurat.integrated)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# Save
save(cds, file=outfile)