library(Seurat)
library(fgsea)
library(dplyr)
library(Matrix)
library(GSEABase)

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
infile5 <- commandArgs(trailingOnly=TRUE)[5]
outfile <- commandArgs(trailingOnly=TRUE)[6]
# infile1 <- 'output/hpbase/DAPT/seurat_annotated_landscaper.RData'
# infile2 <- 'plot/hpbase/DAPT/Landscaper/Allstates_major_group.tsv'
# infile3 <- 'plot/hpbase/DAPT/Landscaper/BIN_DATA'
# infile4 <- 'plot/hpbase/DAPT/Landscaper/Basin.tsv'
# infile5 <- 'data/go_bp_hpbase.RData'

# Gene sets -> fgsea pathways
to_fgsea_pathways <- function(gmt_go) {
  gs_list <- if (inherits(gmt_go, "GeneSetCollection")) as.list(gmt_go) else gmt_go
  pathways <- lapply(gs_list, geneIds)
  nm <- vapply(gs_list, setName, character(1))
  nm[nm == ""] <- vapply(gs_list, setIdentifier, character(1))[nm == ""]
  names(pathways) <- nm
  pathways
}

# Load Seurat
load(infile1)
md <- seurat.integrated@meta.data

# Load GO
objname <- load(infile5)
gmt_go <- get(objname)
rm(list = objname)
gos <- to_fgsea_pathways(gmt_go)

# Read BIN_DATA
bin <- as.matrix(read.table(infile3, header=FALSE, sep="", stringsAsFactors=FALSE))
storage.mode(bin) <- "integer"
rownames(bin) <- rownames(md)
K <- ncol(bin)

# Read Allstates_major_group (pattern + state_id)
allst_df <- read.table(infile2, header=FALSE, sep="", stringsAsFactors=FALSE)
pat_all <- as.matrix(allst_df[, 1:K])
storage.mode(pat_all) <- "integer"
state_id <- allst_df[, K+1]

# Map each cell -> state_id
key_all <- apply(pat_all, 1, paste, collapse=",")
key_bin <- apply(bin, 1, paste, collapse=",")
idx <- match(key_bin, key_all)
cell_state_id <- state_id[idx]

# Basin ids (複数 basin を想定：列が複数なら全部読む)
basin_ids <- as.matrix(read.table(infile4, header=FALSE))
storage.mode(basin_ids) <- "integer"

# Expression matrix
expr <- GetAssayData(seurat.integrated, assay="RNA", layer="data")

# ---- BasinごとにGSEA ----
res_list <- vector("list", nrow(basin_ids))

for (b in seq_len(nrow(basin_ids))) {
  ids <- basin_ids[b,]
  ids <- ids[is.finite(ids)]
  ids <- ids[ids > 0]

  is_basin <- cell_state_id %in% ids

  mu_basin    <- Matrix::rowMeans(expr[, is_basin, drop=FALSE])
  mu_nonbasin <- Matrix::rowMeans(expr[, !is_basin, drop=FALSE])

  score <- mu_basin - mu_nonbasin
  score <- score[is.finite(score)]
  score <- sort(score, decreasing = TRUE)

  set.seed(1)
  fg <- fgsea(
    pathways = gos,
    stats = score,
    minSize = 5,
    maxSize = 500,
    nperm = 10000
  ) %>% arrange(padj, desc(abs(NES)))

  res_list[[b]] <- list(
    basin_ids = ids,
    n_basin = sum(is_basin),
    n_nonbasin = sum(!is_basin),
    score = score,
    fg = fg
  )
}

names(res_list) <- paste0("basin_", seq_along(res_list))

save(res_list, gos, file = outfile)
