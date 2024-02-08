source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
xlsx_sheet <- read_excel(infile2, sheet=1)

# rm NA.
xlsx_sheet <- xlsx_sheet[, 1:7]

# Cluster ID => Celltype name
l <- length(table(Idents(seurat.integrated)))
cols <- c()
for(i in seq(l)){
    clusterid <- as.character(i-1)
    germlayer <- xlsx_sheet$ecto_endo_meso[i]
    cmd1 <- paste0(
        "seurat.integrated <- RenameIdents(seurat.integrated, ",
        "'", clusterid, "' = ", "'", germlayer, "')")
    eval(parse(text=cmd1))
}
Idents(seurat.integrated) <- factor(Idents(seurat.integrated),
    levels=c("Ectoderm", "Mesoderm", "Endoderm", "NA"))
cols <- germlayer_colors

# Save
save(seurat.integrated, file=outfile)
