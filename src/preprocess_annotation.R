source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
annotation <- read_excel(infile, sheet=1)
annotation <- as.data.frame(annotation)

# Save
save(annotation, file=outfile)