source("src/Functions.R")

# Parameter
infile <- 'plot/hpbase/integrated_cov/sce.RData'
outfile <- 'plot/hpbase/cont_DAPT_cov/energy_diff.png'

# Loading
load(infile)

# ===== 1. 統合SCEオブジェクト（sce）に対しhexbin計算 =====
# sce@metadata$hexbin[[1]]：各細胞のhexbin ID
# sce@metadata$hexbin[[2]]：各hexbinの位置情報（data.frame形式、列名例: hexID, x, y, count）
hex_ids <- sce@metadata$hexbin[[1]]
hex_positions <- as.data.frame(sce@metadata$hexbin[[2]])
hex_positions <- rownames_to_column(hex_positions, var = "hexID")

# ===== 2. 各細胞のhexbin ID, energy, condition情報をまとめ、グループ別中央値を計算 =====
df <- data.frame(
  hexID = hex_ids,
  energy = sce$cont_DAPT_energy,
  condition = sce$condition
)
df$hexID <- as.character(df$hexID)
# 3. 元のグリッド上のユニークなhexbin番号を昇順に並べる
unique_hex <- sort(unique(hex_ids))
# 4. df の hexID を、unique_hex 中の位置に置き換える
df$hexID <- match(hex_ids, unique_hex)  # これで値は1～1356になる

# hexbinごと・conditionごとにエネルギーの中央値を計算し、ワイド形式に変換
median_wide <- df %>% 
  group_by(hexID, condition) %>% 
  summarize(med_energy = median(energy, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = condition, values_from = med_energy) %>%
  filter(!is.na(cont) & !is.na(DAPT)) %>%
  mutate(diff = DAPT - cont)

# ===== 3. 六角形の位置情報と結合 =====
hex_diff <- merge(median_wide, hex_positions, by = "hexID", all = FALSE)

# ===== 4. 各プロットの作成 =====

# プロットA：Group A の中央値プロット
p1 <- ggplot(hex_diff, aes(x = x, y = y, fill = cont)) +
  geom_hex(stat = "identity", na.rm = TRUE) +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Control Energy (Median)", fill = "cont")

# プロットB：Group B の中央値プロット
p2 <- ggplot(hex_diff, aes(x = x, y = y, fill = DAPT)) +
  geom_hex(stat = "identity", na.rm = TRUE) +
  scale_fill_viridis_c() +
  coord_fixed() +
  theme_minimal() +
  labs(title = "DAPT Energy (Median)", fill = "DAPT")

# プロットDiff：差分プロット（B – A）
p3 <- ggplot(hex_diff, aes(x = x, y = y, fill = diff)) +
  geom_hex(stat = "identity", na.rm = TRUE) +
  scale_fill_gradient2(low = "blue", mid = "gray80", high = "red", midpoint = 0) +
  coord_fixed() +
  theme_minimal() +
  labs(title = "Diff (DAPT - cont)", fill = "DAPT - cont")

# ===== 5. 3つのプロットを横並びに表示 =====
png(file=outfile, width=1800, height=600)
grid.arrange(p1, p2, p3, ncol = 3)
dev.off()
