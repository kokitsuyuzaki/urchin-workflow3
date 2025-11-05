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
# infile1 <- 'output/hpbase/integrated/seurat_annotated.RData'
# infile2 <- 'output/hpbase/cont/seurat_annotated.RData'
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

## Only Ectoderm in 24h, 36h, 48h samples
target1 <- which(seurat.integrated@meta.data$germlayer == "Ectoderm")
target2 <- grep("24h|36h|48h", seurat.integrated@meta.data$sample)
seurat.integrated <- seurat.integrated[, intersect(target1, target2)]

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

################################################
################## metropolis ##################
################################################
# none
p_bg <- base_plot(umap, mode = "none", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile1, p, width = 7, height = 6, dpi = 300)

# germlayer
p_bg <- base_plot(umap, mode = "germlayer", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile2, p, width = 7, height = 6, dpi = 300)

# cluster
p_bg <- base_plot(umap, mode = "cluster", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile3, p, width = 7, height = 6, dpi = 300)

# sample
p_bg <- base_plot(umap, mode = "sample", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile4, p, width = 7, height = 6, dpi = 300)

# celltype
p_bg <- base_plot(umap, mode = "celltype", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile5, p, width = 7, height = 6, dpi = 300)

# state
p_bg <- base_plot(umap, mode = "state", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile6, p, width = 7, height = 6, dpi = 300)

# energy
p_bg <- base_plot(umap, mode = "energy", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_m) > 0) geom_segment(data = df_m,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_m$col, alpha = scales::rescale(df_m$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile7, p, width = 7, height = 6, dpi = 300)


################################################
################## glauber ##################
################################################
# none
p_bg <- base_plot(umap, mode = "none", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile8, p, width = 7, height = 6, dpi = 300)

# germlayer
p_bg <- base_plot(umap, mode = "germlayer", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile9, p, width = 7, height = 6, dpi = 300)

# cluster
p_bg <- base_plot(umap, mode = "cluster", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile10, p, width = 7, height = 6, dpi = 300)

# sample
p_bg <- base_plot(umap, mode = "sample", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile11, p, width = 7, height = 6, dpi = 300)

# celltype
p_bg <- base_plot(umap, mode = "celltype", meta = meta, E = E, germlayer_cols = NULL)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile12, p, width = 7, height = 6, dpi = 300)

# state
p_bg <- base_plot(umap, mode = "state", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile13, p, width = 7, height = 6, dpi = 300)

# energy
p_bg <- base_plot(umap, mode = "energy", meta = meta, E = E, germlayer_cols = NULL, state_idx=state_idx)

p <- p_bg +
  { if (nrow(df_g) > 0) geom_segment(data = df_g,
      aes(x = x, y = y, xend = xend, yend = yend),
      arrow = arrow(length = unit(0.15, "cm")),
      color = df_g$col, alpha = scales::rescale(df_g$mass, to = c(0.3, 1)),
      linewidth = 0.6) }

ggsave(outfile14, p, width = 7, height = 6, dpi = 300)
