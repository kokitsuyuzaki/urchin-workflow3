source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
celltype_label <- read_excel(infile2, sheetIndex=1)

# rm NA.
celltype_label <- celltype_label[, 1:6]

# Cluster ID => Celltype name
l <- length(table(Idents(seurat.integrated)))
cols <- c()
for(i in seq(l)){
    clusterid <- as.character(i-1)
    celltype <- celltype_label$celltype[i]
    cmd1 <- paste0(
        "seurat.integrated <- RenameIdents(seurat.integrated, ",
        "'", clusterid, "' = ", "'", celltype, "')")
    cmd2 <- paste0(
        "cols <- c(cols, rgb(", celltype_label$R[i],
            ",", celltype_label$G[i],
            ",", celltype_label$B[i], ", maxColorValue=255))")
    eval(parse(text=cmd1))
    eval(parse(text=cmd2))
}
names(cols) <- celltype_label$celltype[seq(l)]
cols <- cols[unique(celltype_label$celltype[seq(l)])]

# Save
save(seurat.integrated, cols, file=outfile)
