source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]

# Loading
load(infile1)
load(infile2)

# Cell cycle-related Genes
s.genes <- intersect(
    rownames(seurat.integrated),
    unlist(lapply(cc.genes$s.genes, function(x){
        annotation$HPU_gene_name[grep(x, annotation$Uniprot_gene_name)]
    })))
g2m.genes <- intersect(
    rownames(seurat.integrated),
    unlist(lapply(cc.genes$g2m.genes, function(x){
        annotation$HPU_gene_name[grep(x, annotation$Uniprot_gene_name)]
    })))
genes_ridgeplot <- genes_ridgeplot_hpbase

# Cellcycle Scoring
seurat.integrated <- CellCycleScoring(seurat.integrated,
    s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Plot
png(file = outfile1, width=600, height=600)
RidgePlot(seurat.integrated, features = genes_ridgeplot, ncol = 3)
dev.off()

# Plot
png(file=outfile2, width=600, height=600)
DimPlot(seurat.integrated, pt.size=2, label.size=6)
dev.off()

# Plot
png(file=outfile3, width=2000, height=1000)
DimPlot(seurat.integrated, split.by="sample", pt.size=2, label.size=6, ncol=5)
dev.off()
