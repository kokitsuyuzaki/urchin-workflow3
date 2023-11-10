source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Extract Group Name
group <- as.character(Idents(seurat.cont))

# Save
write.table(group, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)