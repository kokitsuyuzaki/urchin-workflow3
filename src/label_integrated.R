source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
xlsx_sheet <- read_excel(infile2, sheet=1)

# Cluster ID => Celltype name
l <- length(table(Idents(seurat.integrated)))

celltype <- rep(0, length=ncol(seurat.integrated))
celltype_colors <- rep(0, length=ncol(seurat.integrated))
for(i in seq(l)){
    target <- which(Idents(seurat.integrated) == as.character(i-1))
    cmd <- paste0("rgb(", xlsx_sheet$R[i], ",",
            xlsx_sheet$G[i], ",",
            xlsx_sheet$B[i], ", maxColorValue=255)")
    celltype[target] <- xlsx_sheet$celltype[i]
    celltype_colors[target] <- eval(parse(text=cmd))
}
seurat.integrated$celltype <- celltype
seurat.integrated$celltype_colors <- celltype_colors

# Cluster ID => Germ Layer
germlayer <- rep(0, length=ncol(seurat.integrated))
for(i in seq(l)){
    target <- which(Idents(seurat.integrated) == as.character(i-1))
    cmd <- paste0("rgb(", xlsx_sheet$R[i], ",",
            xlsx_sheet$G[i], ",",
            xlsx_sheet$B[i], ", maxColorValue=255)")
    germlayer[target] <- xlsx_sheet$ecto_endo_meso[i]
}
seurat.integrated$germlayer <- germlayer
seurat.integrated$germlayer_colors <- germlayer_colors[germlayer]

# Save
save(seurat.integrated, file=outfile)
