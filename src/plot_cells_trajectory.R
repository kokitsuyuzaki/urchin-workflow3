source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Plot
png(file=outfile, width=600, height=600)
plot_cells(cds, label_groups_by_cluster=FALSE,
    label_leaves=FALSE, label_branch_points=FALSE)
dev.off()
