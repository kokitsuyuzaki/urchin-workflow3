source("src/Functions2.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile)

# Stratification
cells_to_keep <- rownames(seurat.integrated@meta.data)[grepl("^DAPT-", seurat.integrated$sample)]
seurat.integrated <- subset(seurat.integrated, cells = cells_to_keep)

# 1. sample列から'cont-'と'DAPT-'を削除
seurat.integrated$cluster_day <- gsub("^DAPT-", "", seurat.integrated$sample)

# data slot を counts に置き換え
DefaultAssay(seurat.integrated) <- "RNA"
seurat.integrated[["RNA"]]@data <- seurat.integrated[["RNA"]]@counts

# germlayerごとに分割して出力
layers <- c("Ectoderm", "Mesoderm", "Endoderm")
outfiles <- c(outfile1, outfile2, outfile3)
for (i in 1:3) {
    # Setting
    layer <- layers[i]
    outfile <- outfiles[i]

    # Stratification
    seurat.sub <- subset(seurat.integrated, germlayer == layer)
    
    # 一時ファイル名の生成
    outfile2 <- gsub(".h5ad$", ".h5Seurat", outfile)
    
    unlink(outfile2)
    SaveH5Seurat(seurat.sub, filename = outfile2)
    Convert(outfile2, dest = "h5ad", overwrite = TRUE)
}
