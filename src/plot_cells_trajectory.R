source("src/Functions.R")
library(ggplot2)

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# Loading
load(infile)

# Preprocess
df <- unique(as.data.frame(colData(cds)[, c("celltype", "celltype_colors")]))

celltype_cols <- as.character(df$celltype_colors)  # "#RRGGBB"
names(celltype_cols) <- as.character(df$celltype)

# 念のため、順序固定
colData(cds)$celltype <- factor(colData(cds)$celltype, levels = names(celltype_cols))

# Plot
png(file=outfile, width=600, height=600)
g <- plot_cells(
  cds,
  color_cells_by = "celltype",
  label_groups_by_cluster = FALSE,
  label_leaves = FALSE,
  label_branch_points = FALSE,
  cell_size = 1.5,
  trajectory_graph_segment_size = 2.5,
  trajectory_graph_color = "black"
) +
  scale_color_manual(values = celltype_cols) +
  theme_void() +
  theme(
    legend.position = "none",
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
print(g)
dev.off()