source("src/Functions.R")

args <- commandArgs(trailingOnly = TRUE)
infile1 <- args[1]
infile2 <- args[2]
infile3 <- args[3]
infile4 <- args[4]
infile5 <- args[5]
infile6 <- args[6]
infile7 <- args[7]
outfile1 <- args[8]
outfile2 <- args[9]
outfile3 <- args[10]
outfile4 <- args[11]
outfile5 <- args[12]
outfile6 <- args[13]
outfile7 <- args[14]
outfile8 <- args[15]
outfile9 <- args[16]
outfile10 <- args[17]
outfile11 <- args[18]
outfile12 <- args[19]
outfile13 <- args[20]
outfile14 <- args[21]
# infile1 <- 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
# infile2 <- 'output/hpbase/cont/seurat_annotated_landscaper.RData'
# infile3 <- 'output/hpbase/cont/sbmfcv/BIN_DATA.tsv'
# infile4 <- 'plot/hpbase/cont/Landscaper/Allstates.tsv'
# infile5 <- 'plot/hpbase/cont/Landscaper/E.tsv'
# infile6 <- 'plot/hpbase/cont/P_metropolis_cg.RData'
# infile7 <- 'plot/hpbase/cont/P_glauber_cg.RData'

# Load
load(infile1)
meta <- seurat.integrated@meta.data
umap <- Embeddings(seurat.integrated, "umap")[, 1:2, drop = FALSE]

load(infile2)

meta <- meta[colnames(seurat.integrated), ]
umap <- umap[colnames(seurat.integrated), ]

BIN <- as.matrix(read.table(infile3, header=FALSE))
Allstates <- as.matrix(read.table(infile4, header=FALSE))
E <- unlist(read.delim(infile5, header = FALSE))
load(infile6)
load(infile7)

# ---- state_index を BIN × Allstates から作る ----
Sdim <- ncol(Allstates)

key <- function(v) paste0(v, collapse = "|")
keys_all  <- apply(Allstates, 1, key)
keys_cell <- apply(BIN,       1, key)
idx_map   <- setNames(seq_len(nrow(Allstates)), keys_all)
state_idx <- unname(idx_map[keys_cell])

df_m <- make_arrow_df(dr_m, mu, umap, arrow_color = "#000000ff")
df_g <- make_arrow_df(dr_g, mu, umap, arrow_color = "#000000ff")

#========================================================
#===================  Metropolis  =======================
#========================================================

# 観測状態ベクトルを格子へ内挿（疎→密）
grid_m <- make_vector_grid(df_m, umap, nx = 30, ny = 30)

layers_m <- if (nrow(grid_m) > 0) list(
  metR::geom_streamline(
    data = grid_m,
    aes(x = x, y = y, dx = dx, dy = dy, alpha = speed),
    L = 1,
    min.L = 1,
    res = 2,
    S = 1.5,
    method = "rk4",  # Runge-Kutta法（より正確）
    arrow = arrow(length = unit(0.2, "cm"), type = "closed", angle = 30),
    lineend = "round",
    linewidth = 0.4,
    n = 13,
    colour = "black"
  ),
  scale_alpha(range = c(0.3, 1), guide = "none")
) else NULL

# 背景
p_bg <- base_plot(umap, mode = "none", meta = meta, E = E)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile1, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "germlayer", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile2, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "cluster", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile3, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "sample", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile4, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "celltype", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile5, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "state", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile6, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "energy", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

# ストリームライン
p <- p_bg + layers_m
ggsave(outfile7, p, width = 7, height = 6, dpi = 300)
#========================================================

#========================================================
#===================  Glauber  =======================
#========================================================

# 観測状態ベクトルを格子へ内挿（疎→密）
grid_g <- make_vector_grid(df_g, umap, nx = 30, ny = 30)

layers_g <- if (nrow(grid_g) > 0) list(
  metR::geom_streamline(
    data = grid_g,
    aes(x = x, y = y, dx = dx, dy = dy, alpha = speed),
    L = 1,
    min.L = 1,
    res = 2,
    S = 1.5,
    method = "rk4",  # Runge-Kutta法（より正確）
    arrow = arrow(length = unit(0.2, "cm"), type = "closed", angle = 30),
    lineend = "round",
    linewidth = 0.4,
    n = 13,
    colour = "black"
  ),
  scale_alpha(range = c(0.3, 1), guide = "none")
) else NULL

# 背景
p_bg <- base_plot(umap, mode = "none", meta = meta, E = E)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile8, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "germlayer", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile9, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "cluster", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile10, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "sample", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile11, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "celltype", meta = meta, E = E, germlayer_cols = NULL)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile12, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "state", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile13, p, width = 7, height = 6, dpi = 300)

# 背景
p_bg <- base_plot(umap, mode = "energy", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

# ストリームライン
p <- p_bg + layers_g
ggsave(outfile14, p, width = 7, height = 6, dpi = 300)
#========================================================
