source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
outfile1 <- commandArgs(trailingOnly=TRUE)[5]
outfile2 <- commandArgs(trailingOnly=TRUE)[6]
# infile1 <- 'plot/hpbase/integrated/Landscaper/E.tsv'
# infile2 <- 'plot/hpbase/integrated/Landscaper/Allstates.tsv'
# infile3 <- 'output/hpbase/integrated/binpca/BIN_DATA.tsv'
# infile4 <- 'output/hpbase/integrated/cov.tsv'

##------------------------------------------------------------------------------
## ロード
##------------------------------------------------------------------------------

E         <- as.numeric(read.table(infile1, header = FALSE)[, 1])
Allstates <- as.matrix(read.table(infile2, header = FALSE))
BIN       <- as.matrix(read.table(infile3, header = FALSE))
cell_time <- as.numeric(read.table(infile4, header = FALSE)[, 1])

P    <- nrow(Allstates)
Sdim <- ncol(Allstates)

if (ncol(BIN) != Sdim) {
  stop("Number of columns in BIN_DATA (", ncol(BIN),
       ") does not match ncol(Allstates) (", Sdim, ").")
}
if (length(cell_time) != nrow(BIN)) {
  stop("Length of cell_time (", length(cell_time),
       ") must equal nrow(BIN_DATA) (", nrow(BIN), ").")
}
if (length(E) != P) {
  stop("Length of E (", length(E), ") and nrow(Allstates) (", P, ") do not match.")
}

##------------------------------------------------------------------------------
## セルごとのパターン → Allstates の行番号 state_idx を求める
##------------------------------------------------------------------------------

key_fun <- function(v) paste0(v, collapse = "|")

keys_all  <- apply(Allstates, 1, key_fun)  # 長さ P
keys_cell <- apply(BIN,       1, key_fun)  # 長さ Ncell

idx_map <- setNames(seq_len(P), keys_all)

state_idx <- unname(idx_map[keys_cell])

if (any(is.na(state_idx))) {
  warning("Some cells could not be matched to any row in Allstates (state_idx has NA). ",
          "パターンがAllstatesに含まれていないセルがあります。")
}

##------------------------------------------------------------------------------
## 状態ごとの平均時間 time_state（長さ P）を作成
##------------------------------------------------------------------------------

time_state <- rep(NA_real_, P)

tmp <- tapply(cell_time, state_idx, mean)  # names(tmp) は state_idx の値（1..P）

idx_states_obs <- as.integer(names(tmp))
time_state[idx_states_obs] <- as.numeric(tmp)

##------------------------------------------------------------------------------
## 転移確率行列の構築（Metropolis / Glauber） with 時系列制約
##------------------------------------------------------------------------------

## 許容幅 tol:
## time_state[j] < time_state[i] - tol を禁止
## 36,48,72,96 なら tol=0 で「厳密に過去は禁止、同時刻は許可」
tol_h <- 0

P_metropolis <- transition_matrix_from_E(
  E          = E,
  S          = Allstates,
  beta       = 1,
  kernel     = "metropolis",
  time_state = time_state,
  tol        = tol_h
)

P_glauber <- transition_matrix_from_E(
  E          = E,
  S          = Allstates,
  beta       = 1,
  kernel     = "glauber",
  time_state = time_state,
  tol        = tol_h
)

##------------------------------------------------------------------------------
## 出力
##------------------------------------------------------------------------------

write.table(
  as.matrix(P_metropolis),
  file      = outfile1,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)

write.table(
  as.matrix(P_glauber),
  file      = outfile2,
  sep       = "\t",
  row.names = FALSE,
  col.names = FALSE,
  quote     = FALSE
)