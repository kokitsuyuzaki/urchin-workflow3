source("src/Functions2.R")

# ---------------- args ----------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript convert_seurat_to_h5ad.R <infile.RData> <outfile.h5ad>")
infile  <- args[1]
outfile <- args[2]

# ---------------- load ----------------
load(infile)  # expects seurat.integrated
if (!exists("seurat.integrated")) stop("Object 'seurat.integrated' not found in infile.")

# ---------------- sanitize reductions/keys (avoid NULL attributes) ----------------
if ("pca" %in% names(seurat.integrated@reductions)) {
  if (is.null(seurat.integrated[["pca"]]@key) || seurat.integrated[["pca"]]@key == "")
    seurat.integrated[["pca"]]@key <- "PC_"
  seurat.integrated[["pca"]]@misc <- list()
}

if ("umap" %in% names(seurat.integrated@reductions)) {
  if (is.null(seurat.integrated[["umap"]]@key) || seurat.integrated[["umap"]]@key == "")
    seurat.integrated[["umap"]]@key <- "UMAP_"
  seurat.integrated[["umap"]]@misc <- list()
}

# clear misc at object-level to avoid NULL attributes
seurat.integrated@misc <- list()

# ---------------- slim the object (only what is needed) ----------------
# keep typical assays if present; adjust as appropriate
keep_assays <- intersect(c("SCT","RNA","integrated"), names(seurat.integrated@assays))
if (length(keep_assays) == 0) keep_assays <- names(seurat.integrated@assays)[1]

seurat.integrated <- DietSeurat(
  seurat.integrated,
  assays     = keep_assays,
  counts     = TRUE,
  data       = TRUE,
  scale.data = FALSE,
  dimreducs  = intersect(c("pca","umap"), names(seurat.integrated@reductions)),
  graphs     = NULL,
  misc       = FALSE
)

# assure DefaultAssay is valid
DefaultAssay(seurat.integrated) <- keep_assays[1]

# ---------------- write h5Seurat then convert to h5ad ----------------
outfile2 <- sub("\\.h5ad$", ".h5Seurat", outfile, ignore.case = TRUE)

if (file.exists(outfile2)) unlink(outfile2)
if (file.exists(outfile))  unlink(outfile)

SaveH5Seurat(seurat.integrated, filename = outfile2, overwrite = TRUE)

# Choose a reasonable assay for export (prefer SCT, else RNA, else first)
export_assay <- intersect(c("SCT","RNA"), keep_assays)
export_assay <- if (length(export_assay)) export_assay[1] else keep_assays[1]

Convert(outfile2, dest = "h5ad", assay = export_assay, overwrite = TRUE)

cat(sprintf("Wrote: %s (assay=%s)\n", outfile, export_assay))