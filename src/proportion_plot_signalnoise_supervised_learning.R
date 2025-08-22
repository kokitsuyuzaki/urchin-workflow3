# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# パッケージ
library(ggplot2)
library(dplyr)
library(scales)

# Loading
load(infile)

wanted   <- c("celltype", "germlayer", "Phase", "seurat_clusters", "time")
unwanted <- c("condition", "dbl.dens", "nCount_SCT", "nFeature_SCT", "percent.mt", "percent.rb")

rename_map <- c(
  "dbl.dens"     = "Doublet",
  "nCount_SCT"   = "No. Counts",
  "nFeature_SCT" = "No. Detected Genes",
  "percent.mt"   = "Mitochondria (%)",
  "percent.rb"   = "Ribosome (%)"
)

cap1 <- function(s) ifelse(nchar(s) > 0, paste0(toupper(substr(s,1,1)), substr(s,2,nchar(s))), s)

# ===== 前処理 =====
df <- results %>%
  transmute(
    name  = as.character(target),
    value = as.numeric(common_skill)
  ) %>%
  mutate(
    value = ifelse(is.finite(value), value, 0),   # NaN/Inf/NA -> 0
    group = case_when(
      name %in% wanted   ~ "Wanted",
      name %in% unwanted ~ "Unwanted",
      TRUE ~ "Unwanted"
    ),
    label = dplyr::recode(name, !!!rename_map, .default = name),
    label = cap1(label)
  ) %>%
  arrange(desc(value))
  
df$label <- factor(df$label, levels = df$label)
df$group <- factor(df$group, levels = c("Wanted","Unwanted"))

# ===== プロット =====
p <- ggplot(df, aes(x = label, y = value, fill = group)) +
  geom_col(width = 0.8) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 0.1)) +
  scale_fill_manual(
    name   = "Category",
    values = c(Wanted = "#2ca02c", Unwanted = "#d62728"),
    breaks = c("Wanted", "Unwanted"),
    labels = c("Wanted", "Unwanted")
  ) +
  labs(x = NULL, y = "Common skill") +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

# Plot
ggsave(outfile, p, width = 7, height = 4.5, dpi = 300)