source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
sample_name <- gsub("output/hpbase/", "", gsub("/seurat.RData", "", infile1))
position <- which(sample_name == sample_names)
celltype_label <- read_excel(infile2, sheet=group_names[position])

# rm NA.
celltype_label <- celltype_label[, 1:7]

# Cluster ID => Celltype name
l <- length(table(Idents(seurat.obj)))
for(i in seq(l)){
    clusterid <- as.character(i-1)
    germlayer <- celltype_label$ecto_endo_meso[i]
    cmd1 <- paste0(
        "seurat.obj <- RenameIdents(seurat.obj, ",
        "'", clusterid, "' = ", "'", germlayer, "')")
    eval(parse(text=cmd1))
}
Idents(seurat.obj) <- factor(Idents(seurat.obj),
    levels=c("Ectoderm", "Mesoderm", "Endoderm", "NA"))
cols <- germlayer_colors

# Output
save(seurat.obj, cols, file=outfile)
