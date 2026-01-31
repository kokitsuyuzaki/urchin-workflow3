source("src/Functions.R")

# Parameter
sample_name <- commandArgs(trailingOnly=TRUE)[1]
infile <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# Sub-sampling
seurat.integrated <- subset(seurat.integrated, subset = sample == sample_name)

# Save
save(seurat.integrated, file=outfile)
