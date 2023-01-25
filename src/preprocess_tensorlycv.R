source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Tensor construction
seurat.integrated@meta.data[, c("seurat_clusters", "sample")] %>%
    .countCategory %>% pivot_wider(., names_from="sample", values_from="no_cells") %>% as.data.frame -> data
rownames(data) <- data[,1]
data <- data[,2:ncol(data)]

# Sort
urchin_tensor <- data[seurat_clusters,,]

# Int => Float
urchin_tensor <- as.matrix(urchin_tensor)
urchin_tensor[] <- as.double(urchin_tensor)

# Save Python binary file
np <- import("numpy")
np$save(outfile, r_to_py(urchin_tensor))
