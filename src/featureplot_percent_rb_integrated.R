source("src/Functions.R")

# Parameter
db <- commandArgs(trailingOnly=TRUE)[1]
infile1 <- commandArgs(trailingOnly=TRUE)[2]
infile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile1 <- commandArgs(trailingOnly=TRUE)[4]
outfile2 <- commandArgs(trailingOnly=TRUE)[5]

# Loading
load(infile1)
load(infile2)

# Extract Ribosome Genes
if(db == "hpbase"){
    target.genes <- intersect(rownames(seurat.integrated),
        annotation$HPU_gene_name[grep("ribosome", annotation$NR_genename)])
}else{
target.genes <- intersect(rownames(seurat.integrated),
    annotation$SPU_gene_name[grep("ribosome", annotation$NR_genename)])
}

# Calculate the Percentage
percent.rb <- colSums(seurat.integrated[target.genes, ]) / colSums(seurat.integrated) * 100
seurat.integrated[["percent.rb"]] <- percent.rb

# Plot
png(file=outfile1, width=600, height=600)
FeaturePlot(seurat.integrated, features="percent.rb",
    reduction = "umap", pt.size=2, label.size=6) + xlim(c(-15,15)) + ylim(c(-15,15))
dev.off()

seuratList <- .stratifySeurat(seurat.integrated, group_names)
png(file=outfile2, width=2000, height=1000)
.panelPlot(seuratList, group_names, "percent.rb")
dev.off()
