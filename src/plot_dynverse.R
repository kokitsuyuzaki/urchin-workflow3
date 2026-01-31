source("src/Functions3.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
outfile3 <- commandArgs(trailingOnly=TRUE)[5]
# infile1 = 'output/hpbase/integrated/seurat_annotated_landscaper.RData'
# infile2 = 'output/hpbase/integrated/dynverse/paga_celltype.RData'

# Load
load(infile1)
load(infile2)

# 追加: cell_id を安全に取り出すヘルパー
get_cell_ids <- function(model) {
  if (!is.null(model$cell_ids)) return(as.character(model$cell_ids))
  if (!is.null(model$dimred) && !is.null(rownames(model$dimred))) return(rownames(model$dimred))
  if (!is.null(model$counts) && !is.null(colnames(model$counts))) return(colnames(model$counts))
  stop("cell ids not found in model (tried $cell_ids, rownames($dimred), colnames($counts)).")
}

# model 側の cell_id
cell_ids_model <- get_cell_ids(model)

# 2) Seurat メタデータから対応表を作る（cell_id = 行名）
meta <- seurat.integrated@meta.data %>%
  tibble::rownames_to_column("cell_id") %>%
  mutate(
    seurat_clusters = as.character(seurat_clusters),
    sample          = as.character(sample),
    celltype        = as.character(celltype)
  )

# 3) model に存在する cell のみに絞り、model 順に並べる
meta_m <- meta %>%
  semi_join(tibble(cell_id = cell_ids_model), by = "cell_id") %>%
  arrange(match(cell_id, cell_ids_model))

# チェック（不足があればここでわかる）
missing_n <- length(setdiff(cell_ids_model, meta$cell_id))
if (missing_n > 0) {
  message("WARNING: model にあるが meta に無い cell が ", missing_n, " 個あります。")
}

# 4) plot_dimred() が受け取れる形（data.frame: cell_id, group_id）を作る
grp_clusters <- meta_m %>% select(cell_id, group_id = seurat_clusters)
grp_sample   <- meta_m %>% select(cell_id, group_id = sample)
grp_celltype <- meta_m %>% select(cell_id, group_id = celltype)

# Plot
print("outfile1")
png(file = outfile1, width = 1200, height = 1200)
g1 <- plot_dimred(model, grouping = grp_clusters)
if (length(g1$layers) > 0) {
  for (i in seq_along(g1$layers)) {
    # Trajectory線を太くする
    if (inherits(g1$layers[[i]]$geom, "GeomPath") || 
        inherits(g1$layers[[i]]$geom, "GeomSegment")) {
      g1$layers[[i]]$aes_params$size <- 5  # 線の太さ
    }
    # データ点を3.0倍に
    if (inherits(g1$layers[[i]]$geom, "GeomPoint")) {
      if (!is.null(g1$layers[[i]]$aes_params$size)) {
        g1$layers[[i]]$aes_params$size <- g1$layers[[i]]$aes_params$size * 3.0
      } else {
        g1$layers[[i]]$aes_params$size <- 3.0  # デフォルトサイズの3.0倍
      }
    }
  }
}
print(g1)
dev.off()

print("outfile2")
png(file = outfile2, width = 1200, height = 1200)
g2 <- plot_dimred(model, grouping = grp_sample)
if (length(g2$layers) > 0) {
  for (i in seq_along(g2$layers)) {
    # Trajectory線を太くする
    if (inherits(g2$layers[[i]]$geom, "GeomPath") || 
        inherits(g2$layers[[i]]$geom, "GeomSegment")) {
      g2$layers[[i]]$aes_params$size <- 5
    }
    # データ点を3.0倍に
    if (inherits(g2$layers[[i]]$geom, "GeomPoint")) {
      if (!is.null(g2$layers[[i]]$aes_params$size)) {
        g2$layers[[i]]$aes_params$size <- g2$layers[[i]]$aes_params$size * 3.0
      } else {
        g2$layers[[i]]$aes_params$size <- 3.0
      }
    }
  }
}
print(g2)
dev.off()

# celltype用の色マッピングを作成
celltype_color_mapping <- seurat.integrated@meta.data %>%
  select(celltype, celltype_colors) %>%
  distinct() %>%
  deframe()

print("outfile3")
png(file = outfile3, width = 1200, height = 1200)

# plot_dimredを実行してデータとtrajectory情報を取得
g3_original <- plot_dimred(model, grouping = grp_celltype)

# データを取得
plot_data <- g3_original$data

# celltypeの色マッピングを適用
plot_data$celltype <- plot_data$color  # 現在のcolor列はcelltype名
plot_data$color_code <- celltype_color_mapping[plot_data$celltype]

# 新しくプロットを作成
library(ggplot2)
g3 <- ggplot(plot_data, aes(x = comp_1, y = comp_2)) +
  geom_point(aes(color = celltype), size = 3) +
  scale_color_manual(values = celltype_color_mapping) +
  coord_equal() +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.major = element_blank(),  # 主グリッド線を削除
    panel.grid.minor = element_blank(),  # 副グリッド線を削除
    panel.border = element_blank(),      # 枠線を削除
    axis.line = element_blank(),          # 軸線を削除
    axis.text.x = element_blank(),        # x軸の数値を削除
    axis.text.y = element_blank(),        # y軸の数値を削除
    axis.ticks = element_blank(),         # 軸の目盛りを削除
    axis.title.x = element_blank(),       # x軸のラベルを削除
    axis.title.y = element_blank()        # y軸のラベルを削除
  )

# trajectory線を元のプロットから抽出して追加（もし必要なら）
if(length(g3_original$layers) > 0) {
  for(layer in g3_original$layers) {
    if(inherits(layer$geom, "GeomPath") || inherits(layer$geom, "GeomSegment")) {
      g3 <- g3 + layer
    }
  }
}

print(g3)
dev.off()