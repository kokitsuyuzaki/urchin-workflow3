source("src/Functions.R")

# Parameter
set.seed(1234)
sample <- commandArgs(trailingOnly=TRUE)[1]
infile1 <- commandArgs(trailingOnly=TRUE)[2]
infile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile <- commandArgs(trailingOnly=TRUE)[4]
indir <- paste0("output/hpbase/", sample, "/outs/filtered_feature_bc_matrix")

# Loading
geneid_to_genename <- read.csv(infile1, row.names=1)
load(infile2)
data <- Read10X(data.dir=indir)

# Gene ID => Gene Name
tmp <- convertRowID(
    input = as.matrix(data),
    rowID = rownames(data),
    LtoR = geneid_to_genename[, c(1,4)])
data <- tmp$output
data <- data[which(rownames(data) != ""), ]
data <- as(data, "sparseMatrix")

# Seurat Object
seurat.obj2 <- CreateSeuratObject(counts = data)

# SCTransform
seurat.obj2 <- SCTransform(seurat.obj2)

# Inheritance from the original object
seurat.obj2[["pca"]] <- seurat.obj[["pca"]]
seurat.obj2[["umap"]] <- seurat.obj[["umap"]]

object.names <- c("orig.ident", "nCount_SCT", "seurat_clusters",
    "celltype_colors", "nCount_RNA", "nFeature_SCT", "sample",
    "germlayer", "nFeature_RNA", "SCT_snn_res.0.8", "celltype",
    "germlayer_colors")
for(i in object.names){
    cmd1 <- paste0("tmp <- seurat.obj$", i)
    cmd2 <- paste0("names(tmp) <- colnames(seurat.obj2)")
    cmd3 <- paste0("seurat.obj2$", i, " <- tmp")
    eval(parse(text=cmd1))
    eval(parse(text=cmd2))
    eval(parse(text=cmd3))
}
Idents(seurat.obj2) <- Idents(seurat.obj)

# Save
seurat.obj <- seurat.obj2
save(seurat.obj, file=outfile)
