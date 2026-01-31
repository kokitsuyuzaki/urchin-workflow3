source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
annotation <- read_excel(infile, sheet=1)

# Preprocesssing
x <- annotation$NR_genename
m <- regexpr("\\(LOC.*\\)", x)
SPU_gene_id <- substr(x, m+1, m + attr(m, "match.length") - 2)

out <- data.frame(
    HPU_gene_id = annotation$HPU_gene_id,
    HPU_gene_name = annotation$HPU_gene_name,
    SPU_gene_name = annotation$SPU_gene_name,
    SPU_gene_id = SPU_gene_id)

# AN-2
out <- rbind(out, c("HPU_NEWLOC1", "Hp-AN-2", "none", "LOC763123"))

# Save
write.csv(out, file=outfile)