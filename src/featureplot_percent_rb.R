source("src/Functions.R")

# Parameter
db <- commandArgs(trailingOnly=TRUE)[1]
infile1 <- commandArgs(trailingOnly=TRUE)[2]
infile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile1)
load(infile2)

# Extract Ribosome Genes
if(db == "hpbase"){
    target.genes <- intersect(rownames(seurat.obj),
        annotation$HPU_gene_name[grep("ribosome", annotation$NR_genename)])
}else{
target.genes <- intersect(rownames(seurat.obj),
    annotation$SPU_gene_name[grep("ribosome", annotation$NR_genename)])
}

# Calculate the Percentage
percent.rb <- colSums(seurat.obj[target.genes, ]) / colSums(seurat.obj) * 100
seurat.obj[["percent.rb"]] <- percent.rb

# Plot
png(file=outfile, width=600, height=600)
FeaturePlot(seurat.obj, features="percent.rb", pt.size=2, label.size=6)
dev.off()
