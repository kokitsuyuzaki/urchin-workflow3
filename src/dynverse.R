source("src/Functions3.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[2]
infile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile <- commandArgs(trailingOnly=TRUE)[4]
meth <- commandArgs(trailingOnly=TRUE)[5]
# infile1 = "output/hpbase/DAPT/seurat_annotated.RData"
# infile2 = "output/hpbase/DAPT/dynverse/start_id_celltype.RData"

# Loading
load(infile1)
load(infile2)

# Preprocess
mat1 <- seurat.integrated[["RNA"]]@counts
mat2 <- seurat.integrated[["SCT"]]@data
mat2@x <- log2(mat2@x + 1)
umap_coords <- seurat.integrated@reductions$umap@cell.embeddings

# Filtering
common.row.names <- intersect(rownames(mat1), rownames(mat2))
mat1 <- t(mat1[common.row.names, ])
mat2 <- t(mat2[common.row.names, ])

# Dynverse Setting
dataset <- wrap_expression(
  counts = mat1,
  expression = mat2
)

dataset <- add_dimred(
    dataset,
    umap_coords)

dataset <- add_prior_information(
  dataset,
  start_id = start_id)

# Trajectory Inference
model <- infer_trajectory(dataset, 
                         method = meth,
                         parameters = list(
                           n_neighbors = 30,      # より多くの近傍を考慮（デフォルト: 15）
                           n_comps = 50,         # より多くのPCA成分を使用
                           n_dcs = 15,           # diffusion components
                           resolution = 1,       # クラスタリング解像度
                           embedding_type = "fa" # or "umap"
                         ))

# Save
save(model, file=outfile)
