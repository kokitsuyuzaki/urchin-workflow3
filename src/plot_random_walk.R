source("src/Functions.R")

# Parameter
args <- commandArgs(trailingOnly=TRUE)
infile1  <- args[1]  # major_group.tsv
infile2  <- args[2]  # random_walk list output (.RData)
infile3  <- args[3]  # Basin.tsv

# outfiles: 10枚（5 celltype × basin/non）
outfile_ab_basin   <- args[4]
outfile_or_basin   <- args[5]
outfile_ci_basin   <- args[6]
outfile_an_basin   <- args[7]
outfile_ne_basin   <- args[8]
outfile_ab_non     <- args[9]
outfile_or_non     <- args[10]
outfile_ci_non     <- args[11]
outfile_an_non     <- args[12]
outfile_ne_non     <- args[13]

# ----------------------------
# Helpers
# ----------------------------
save_blank_png <- function(outfile, width = 8, height = 6, dpi = 300) {
  png(outfile, width = width, height = height, units = "in", res = dpi)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  dev.off()
}

# ----------------------------
# Load labels (states)
# ----------------------------
major_group <- read.table(
  infile1, header = FALSE, sep = "", quote = "",
  stringsAsFactors = FALSE, check.names = FALSE
)

# last column is label
label_col <- ncol(major_group)
labels_raw <- as.character(major_group[[label_col]])
# drop suffix like ".1", ".2" for base celltype matching
labels_base <- trimws(sub("\\..*$", "", labels_raw))

n_states <- length(labels_base)
row_id <- seq_len(n_states)

# ----------------------------
# Basin: interpret as ROW INDICES
# ----------------------------
basin_rows <- scan(infile3, what = integer(), quiet = TRUE)
basin_rows <- basin_rows[!is.na(basin_rows)]
basin_rows <- unique(basin_rows)
# keep only valid indices
basin_rows <- basin_rows[basin_rows >= 1 & basin_rows <= n_states]

is_basin_row <- row_id %in% basin_rows

# ----------------------------
# Load random-walk probability list
# ----------------------------
load(infile2)
if (exists("P_g_list")) {
  P_list <- P_g_list
} else if (exists("P_m_list")) {
  P_list <- P_m_list
} else {
  stop("P_g_list または P_m_list が見つかりません")
}

if (!is.list(P_list) || length(P_list) == 0) stop("P_list is empty or not a list")

# Steps (0..n)
steps <- 0:(length(P_list) - 1)

# Destination grouping
dest_levels <- c("Aboral_ectoderm", "Oral_ectoderm", "Ciliary_band", "Anus", "Neurons", "Other")
focus_types <- c("Aboral_ectoderm", "Oral_ectoderm", "Ciliary_band", "Anus", "Neurons")
celltype2 <- ifelse(labels_base %in% focus_types, labels_base, "Other")

# Colors (destination colors)
pal <- c(
  "Aboral_ectoderm" = "#008080",
  "Oral_ectoderm"   = "#FFFF00",
  "Ciliary_band"    = "#00FF26",
  "Anus"            = "#3F3F7F",
  "Neurons"         = "#FF00FF",
  "Other"           = "grey70"
)

plot_one_start <- function(start_label, basin_flag, outfile) {
  idx_start <- which(labels_base == start_label & is_basin_row == basin_flag)

  # If start rows not found -> blank
  if (length(idx_start) == 0) {
    save_blank_png(outfile)
    return(invisible(NULL))
  }

  df_list <- vector("list", length(P_list))

  for (k in seq_along(P_list)) {
    Pk <- P_list[[k]]
    if (!is.matrix(Pk)) stop(sprintf("P_list[[%d]] is not a matrix", k))

    # subset start rows
    P_dest_rows <- Pk[idx_start, , drop = FALSE]
    p_dest <- colMeans(P_dest_rows)

    # aggregate probabilities by destination label group
    tmp <- tapply(p_dest, celltype2, sum)

    v <- setNames(rep(0, length(dest_levels)), dest_levels)
    hit <- intersect(names(tmp), dest_levels)
    if (length(hit) > 0) v[hit] <- as.numeric(tmp[hit])

    df_list[[k]] <- tibble(
      step = steps[k],
      Aboral_ectoderm = v["Aboral_ectoderm"],
      Oral_ectoderm   = v["Oral_ectoderm"],
      Ciliary_band    = v["Ciliary_band"],
      Anus            = v["Anus"],
      Neurons         = v["Neurons"],
      Other           = v["Other"]
    )
  }

  df <- bind_rows(df_list) %>%
    pivot_longer(-step, names_to = "dest", values_to = "prob") %>%
    mutate(dest = factor(dest, levels = dest_levels))

  # If everything is zero -> blank
  if (all(df$prob == 0)) {
    save_blank_png(outfile)
    return(invisible(NULL))
  }

  p <- ggplot(df, aes(x = step, y = prob, color = dest)) +
    geom_line(linewidth = 2.0) +
    geom_point(size = 2.5) +
    scale_color_manual(values = pal, guide = "none") +
    scale_x_continuous(breaks = steps) +
    coord_cartesian(ylim = c(0, 1)) +
    labs(x = "Step", y = "Probability") +
    theme_classic(base_size = 12) +
    theme(text = element_text(size = 36))

  ggsave(outfile, p, width = 8, height = 6, dpi = 300)
}

# ----------------------------
# Run (5 starts × basin/nonbasin)
# ----------------------------
plot_one_start("Aboral_ectoderm", TRUE,  outfile_ab_basin)
plot_one_start("Oral_ectoderm",   TRUE,  outfile_or_basin)
plot_one_start("Ciliary_band",    TRUE,  outfile_ci_basin)
plot_one_start("Anus",            TRUE,  outfile_an_basin)
plot_one_start("Neurons",         TRUE,  outfile_ne_basin)

plot_one_start("Aboral_ectoderm", FALSE, outfile_ab_non)
plot_one_start("Oral_ectoderm",   FALSE, outfile_or_non)
plot_one_start("Ciliary_band",    FALSE, outfile_ci_non)
plot_one_start("Anus",            FALSE, outfile_an_non)
plot_one_start("Neurons",         FALSE, outfile_ne_non)