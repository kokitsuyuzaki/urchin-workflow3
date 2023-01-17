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

# Cell cycle-related Genes
if(db == "hpbase"){
    s.genes <- intersect(
        rownames(seurat.obj),
        unlist(lapply(cc.genes$s.genes, function(x){
            annotation$HPU_gene_name[grep(x, annotation$Uniprot_gene_name)]
        })))
    g2m.genes <- intersect(
        rownames(seurat.obj),
        unlist(lapply(cc.genes$g2m.genes, function(x){
            annotation$HPU_gene_name[grep(x, annotation$Uniprot_gene_name)]
        })))
    genes_ridgeplot <- genes_ridgeplot_hpbase
}else{
    s.genes <- intersect(
        rownames(seurat.obj),
        unlist(lapply(cc.genes$s.genes, function(x){
            annotation$SPU_gene_name[grep(x, annotation$Uniprot_gene_name)]
        })))
    g2m.genes <- intersect(
        rownames(seurat.obj),
        unlist(lapply(cc.genes$g2m.genes, function(x){
            annotation$SPU_gene_name[grep(x, annotation$Uniprot_gene_name)]
        })))
    s.genes <- setdiff(s.genes, "none")
    g2m.genes <- setdiff(g2m.genes, "none")
    genes_ridgeplot <- genes_ridgeplot_echinobase
}

# Cellcycle Scoring
seurat.obj <- try(CellCycleScoring(seurat.obj,
    s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE))

if("try-error" %in% is(seurat.obj)){
    file.create(outfile1)
    file.create(outfile2)
    q("no")
}

# Plot
png(file=outfile1, width=600, height=600)
RidgePlot(seurat.obj, features=genes_ridgeplot, ncol=3)
dev.off()
# Plot
png(file=outfile2, width=600, height=600)
DimPlot(seurat.obj, pt.size=2, label.size=6)
dev.off()
