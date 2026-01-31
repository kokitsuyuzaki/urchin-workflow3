source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly = TRUE)[1]
infiles <- paste0(
  "output/hpbase/an2/",
  c("cont-24h", "cont-36h", "cont-48h", "cont-72h", "cont-96h",
    "DAPT-24h", "DAPT-36h", "DAPT-48h", "DAPT-72h", "DAPT-96h"),
  "/seurat.RData"
)

# -------------------------
# 1) load per-sample objects
# -------------------------
seurat.list <- lapply(infiles, function(x) {
  load(x)  # seurat.obj
  seurat.obj
})

# integrated object
load("output/hpbase/integrated/seurat_annotated.RData")  # seurat.integrated
stopifnot(exists("seurat.integrated"))

# -------------------------
# 2) collect Hp-AN-2 expression and align to integrated cells
# -------------------------
gene <- "Hp-AN-2"

an2.list <- lapply(seq_along(seurat.list), function(i) {
  obj <- seurat.list[[i]]
  expr <- FetchData(obj, vars = gene)[, 1]
  names(expr) <- paste0(colnames(obj), "_", i)
  expr
})
an2.vec <- unlist(an2.list)

cells <- colnames(seurat.integrated)
an2.aligned <- an2.vec[cells]
if (anyNA(an2.aligned)) stop("Hp-AN-2 alignment failed: NA found. Check cell naming.")

# -------------------------
# 3) get RNA counts (v4/v5対応) and add Hp-AN-2 as a gene row
# -------------------------
DefaultAssay(seurat.integrated) <- "RNA"

rna <- seurat.integrated[["RNA"]]

# ===== v4 (Assay): そのまま counts を取る =====
counts_all <- GetAssayData(seurat.integrated, assay = "RNA", slot = "counts")
# 念のため
stopifnot(all(colnames(counts_all) == cells))

# ----- Hp-AN-2 を counts_all に追加 -----
if (gene %in% rownames(counts_all)) {
  counts_all[gene, ] <- as.numeric(an2.aligned)
} else {
  new_row <- Matrix(
    as.numeric(an2.aligned),
    nrow = 1,
    dimnames = list(gene, colnames(counts_all))
  )
  counts_all <- rbind(counts_all, new_row)
}

# -------------------------
# 4) replace RNA assay with a clean single-layer assay (v4形式) containing Hp-AN-2
# -------------------------
new_rna <- CreateAssayObject(counts = counts_all)
seurat.integrated[["RNA"]] <- new_rna
DefaultAssay(seurat.integrated) <- "RNA"

# -------------------------
# 5) Save Seurat object
# -------------------------
save(seurat.integrated, file = outfile)
