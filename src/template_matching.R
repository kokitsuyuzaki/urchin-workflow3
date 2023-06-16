source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
template <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile)

# Template Matching
all_genes <- as.matrix(seurat.integrated@assays$RNA@counts)
template_gene <- all_genes[template, ]
other_genes <- all_genes[
    setdiff(seq(nrow(all_genes)),
        which(rownames(all_genes) == template)), ]

cor_other_genes <- t(apply(other_genes, 1, function(x){
    out <- cor.test(x, template_gene)
    c(out$estimate, out$p.value)
}))

bh <- p.adjust(cor_other_genes[,2], method="BH")
qval <- qvalue(cor_other_genes[,2])

cor_other_genes <- cbind(rownames(cor_other_genes),
    cor_other_genes, bh, qval$lfdr, qval$qvalues)
colnames(cor_other_genes) <- c("Symbol", "Cor", "p-value", "FDR_BH", "FDR_LFDR", "FDR_Qvalue")

# Output
save(cor_other_genes, file=outfile1)
write_xlsx(as.data.frame(cor_other_genes), outfile2)
