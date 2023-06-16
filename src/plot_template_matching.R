source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]
template <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile)

# Plot
png(file=outfile, width=1000, height=700)
layout(rbind(1:3, 4:6))
hist(as.numeric(cor_other_genes[,2]),
    main="Pearson Correlation Coefficient",
    breaks=100)
hist(as.numeric(cor_other_genes[,3]),
    main="P-value",
    breaks=100)
hist(as.numeric(cor_other_genes[,4]),
    main="FDR (BH)",
    breaks=100)
hist(as.numeric(cor_other_genes[,5]),
    main="FDR (LFDR)",
    breaks=100)
hist(as.numeric(cor_other_genes[,6]),
    main="FDR (Q-value)",
    breaks=100)
dev.off()
