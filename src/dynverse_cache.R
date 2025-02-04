source("src/Functions3.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]
meth <- commandArgs(trailingOnly=TRUE)[2]

# Singularity Image Fileの取得
eval(parse(text=paste0("ti_", meth, "()")))

# Save
file.create(outfile)
