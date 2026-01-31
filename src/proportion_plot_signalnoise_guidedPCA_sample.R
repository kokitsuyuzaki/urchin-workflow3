# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# パッケージ
library(ggplot2)
library(dplyr)
library(scales)

# Loading
load(infile)  # <- 'results' に contrib_grouped_frac が入っている前提

# 1) ベクトル → データフレーム
vals <- results$contrib_grouped_frac
df <- data.frame(name = names(vals), value = as.numeric(vals), stringsAsFactors = FALSE)

# 2) 区分（Wanted/Unwanted）
wanted   <- c("Phase", "germlayer", "time", "celltype", "seurat_clusters")
unwanted <- c("nCount_SCT", "dbl.dens", "percent.mt", "percent.rb", "nFeature_SCT", "condition")

df <- df %>%
  mutate(
    value = ifelse(is.finite(value), value, 0),   # 安全策
    group = case_when(
      name %in% wanted   ~ "Wanted",
      name %in% unwanted ~ "Unwanted",
      TRUE               ~ "Unwanted"
    )
  )

# 3) 表示名の書き換え＋先頭1文字だけ大文字化
rename_map <- c(
  "dbl.dens"     = "Doublet",
  "nCount_SCT"   = "No. Counts",
  "nFeature_SCT" = "No. Detected Genes",
  "percent.mt"   = "Mitochondria (%)",
  "percent.rb"   = "Ribosome (%)"
)

cap1 <- function(s) ifelse(nchar(s) > 0, paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s))), s)

df <- df %>%
  mutate(label = ifelse(name %in% names(rename_map), rename_map[name], name),
         label = cap1(label))

# ★★★ 特定の項目を明示的に除外 ★★★
df <- df %>% filter(!name %in% c("time", "condition"))

# 4) wanted, unwantedの定義順に並べ替え
order_vec <- c(wanted, unwanted)
order_vec <- order_vec[order_vec %in% df$name]  # データに存在するもののみ
df$name <- factor(df$name, levels = order_vec)
df <- df %>% arrange(name)
# coord_flip()では下から上に並ぶので、逆順にしてlabelをfactor化
df$label <- factor(df$label, levels = rev(df$label))
df$group <- factor(df$group, levels = c("Wanted", "Unwanted"))

# 5) プロット
p <- ggplot(df, aes(x = label, y = value, fill = group)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = function(x) x * 100) +
  scale_fill_manual(
    name   = "",
    values = c(Wanted = "#2ca02c", Unwanted = "#d62728"),
    breaks = c("Wanted", "Unwanted"),
    labels = c("Wanted", "Unwanted")
  ) +
  labs(x = NULL, y = "Explained Variance Fraction (%)") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 16),
    legend.text = element_text(size = 16)
  )

# Plot
ggsave(outfile, p, width = 7, height = 4.5, dpi = 300)
