source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]
# infile1 <- 'plot/hpbase/integrated/Landscaper/major_group.tsv'
# infile2 <- 'plot/hpbase/integrated/P_glauber_rw.tsv'

# Load
major_group <- read.table(infile1, header = FALSE, sep = "", quote = "",
                          stringsAsFactors = FALSE, check.names = FALSE)
P_rw <- as.matrix(read.delim(infile2, header=FALSE, sep="\t"))

## ラベル整形（末尾 .1 など除去＋空白除去）
labels_base <- trimws(sub("\\..*$", "", as.character(major_group[[ncol(major_group)]])))

## 対象（例：Oral_ectoderm 列）  ※ここはあなたの条件に合わせて
idx_OE <- which(labels_base == "Oral_ectoderm")
stopifnot(length(idx_OE) >= 1)
mat_OE  <- P_rw[, idx_OE, drop = FALSE]
mean_OE <- rowMeans(mat_OE)

## 生データ → 不要ラベル除外 → 同名 celltype を集約（平均±SD）
df <- tibble(celltype = labels_base, val = mean_OE) %>%
  filter(!grepl("(?i)unknown|uncharacterized", celltype)) %>%
  group_by(celltype) %>%
  summarise(mean_prob = mean(val),
            sd_prob   = if (n()>1) sd(val) else NA_real_,
            .groups = "drop") %>%
  arrange(desc(mean_prob)) %>%
  mutate(celltype = factor(celltype, levels = celltype))  # 並びを固定

## 色：強調のみ
df$fill_key <- ifelse(df$celltype %in% c("Ciliary_band","Neurons"),
                      as.character(df$celltype), "Other")
pal <- c("Ciliary_band"="#E64B35", "Neurons"="#4DBBD5", "Other"="grey70")

## プロット（降順・タイトルなし）
p <- ggplot(df, aes(x = celltype, y = mean_prob, fill = fill_key)) +
  geom_col(width = 0.8) +
  { if (!all(is.na(df$sd_prob)))
      geom_errorbar(aes(ymin = pmax(0, mean_prob - sd_prob),
                        ymax = mean_prob + sd_prob),
                    width = 0.3, alpha = 0.6) } +
  scale_fill_manual(values = pal, breaks = c("Ciliary_band","Neurons","Other"), name = NULL) +
  scale_x_discrete(limits = levels(df$celltype)) +   # 念のため固定
  labs(x = NULL, y = "Probability") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right")

ggsave(outfile, p, width = 10, height = 6, dpi = 300)