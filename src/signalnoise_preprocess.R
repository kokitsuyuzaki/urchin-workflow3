source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile1)
load(infile2)
load(infile3)

# Time
seurat.integrated[["time"]] <- gsub("DAPT-", "", gsub("cont-", "", seurat.integrated@meta.data$sample))
# Experimental Conditin
seurat.integrated[["condition"]] <- gsub("-.*h", "", seurat.integrated@meta.data$sample)

# Extract Mitochondria Genes
target.genes <- intersect(rownames(seurat.integrated),
    annotation$HPU_gene_name[grep("mitochondria", annotation$NR_genename)])

# Calculate the Percentage
percent.mt <- colSums(seurat.integrated[target.genes, ]) / colSums(seurat.integrated) * 100
seurat.integrated[["percent.mt"]] <- percent.mt

# Extract Ribosome Genes
target.genes <- intersect(rownames(seurat.integrated),
    annotation$HPU_gene_name[grep("ribosome", annotation$NR_genename)])

# Calculate the Percentage
percent.rb <- colSums(seurat.integrated[target.genes, ]) / colSums(seurat.integrated) * 100
seurat.integrated[["percent.rb"]] <- percent.rb

# Cellcycle
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

# Cellcycle Scoring
seurat.integrated <- CellCycleScoring(seurat.integrated,
    s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Doublet Finder
seurat.integrated[["dbl.dens"]] <-dbl.dens

# Save
save(seurat.integrated, file = outfile)


