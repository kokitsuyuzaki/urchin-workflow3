# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile <- commandArgs(trailingOnly=TRUE)[2]

# パッケージ
library(Seurat)
library(ranger)
library(caret)
library(dplyr)
library(tibble)
library(purrr)

# Loading
load(infile)

# ---- 入力: Seuratオブジェクト ----
obj <- seurat.integrated

# ---- 特徴量：SCTの高変動遺伝子（上位2000） ----
assay_use <- "SCT"
if (is.null(VariableFeatures(obj, assay = assay_use)) ||
    length(VariableFeatures(obj, assay = assay_use)) == 0) {
  obj <- FindVariableFeatures(obj, assay = assay_use, nfeatures = 3000)
}
hvg <- head(VariableFeatures(obj, assay = assay_use), 2000)

# セル×遺伝子の行列（denseにする。メモリに注意）
X <- t(as.matrix(GetAssayData(obj, assay = assay_use, slot = "data")[hvg, , drop = FALSE]))

# メタデータ
Y <- obj@meta.data

# ---- 予測対象（質問で挙げた列）----
targets <- c("seurat_clusters", "celltype", "germlayer", "time", "condition",
             "nFeature_SCT", "nCount_SCT", "percent.mt", "percent.rb","Phase", "dbl.dens")
Y <- Y[, targets, drop = FALSE]
Y[, "celltype"] <- factor(Y[, "celltype"])
Y[, "germlayer"] <- factor(Y[, "germlayer"])
Y[, "time"] <- factor(Y[, "time"])
Y[, "condition"] <- factor(Y[, "condition"])
Y[, "Phase"] <- factor(Y[, "Phase"])

.dummify_with_groups <- function(df) {
  if (!is.data.frame(df)) df <- as.data.frame(df)
  n <- NROW(df)

  mats <- list()
  grps <- list()

  add_col <- function(mat, group) {
    mats[[length(mats) + 1]] <<- mat
    grps[[length(grps) + 1]] <<- rep(group, ncol(mat))
  }

  for (nm in names(df)) {
    x <- df[[nm]]

    # ---- 型ごとの前処理（NAは落とさず埋める）----
    if (is.factor(x) || is.character(x)) {
      f <- as.factor(x)  # NAはNAのまま
      lev <- levels(f)
      k <- length(unique(f[!is.na(f)]))  # 非NA水準数

      if (k == 0) {
        next  # 全NAならスキップ（0列）
      } else if (k == 1) {
        # 1水準: 非NAなら1, NAなら0 の定数列
        v <- as.integer(!is.na(f))
        col <- matrix(v, nrow = n, ncol = 1)
        colnames(col) <- paste0(nm, "_", lev[!is.na(lev)][1] %||% "LEVEL")
        add_col(col, nm)
      } else {
        # 2水準以上: 完全手実装のワンホット（NA行は全0）
        lev_nonna <- lev[lev %in% unique(f[!is.na(f)])]
        k2 <- length(lev_nonna)
        M <- matrix(0L, nrow = n, ncol = k2)
        for (j in seq_len(k2)) {
          M[, j] <- as.integer(!is.na(f) & f == lev_nonna[j])
        }
        colnames(M) <- paste0(nm, "_", lev_nonna)
        add_col(M, nm)
      }

    } else if (is.logical(x)) {
      v <- as.integer(replace(x, is.na(x), FALSE))
      col <- matrix(v, nrow = n, ncol = 1)
      colnames(col) <- nm
      add_col(col, nm)

    } else if (inherits(x, "Date") || inherits(x, "POSIXt") || inherits(x, "difftime")) {
      v <- as.numeric(x); v[is.na(v)] <- 0
      if (length(v) != n) stop(sprintf("列'%s'の長さが不正", nm))
      col <- matrix(v, nrow = n, ncol = 1); colnames(col) <- nm
      add_col(col, nm)

    } else {
      # 数値その他
      if (is.matrix(x)) {
        if (!(nrow(x) == n && ncol(x) == 1))
          stop(sprintf("列'%s'が多列/長さ不一致です。1列にしてから渡してください。", nm))
        v <- as.numeric(x[, 1])
      } else {
        v <- as.numeric(x)
      }
      v[is.na(v)] <- 0
      if (length(v) != n) stop(sprintf("列'%s'の長さが不正", nm))
      col <- matrix(v, nrow = n, ncol = 1); colnames(col) <- nm
      add_col(col, nm)
    }
  }

  if (length(mats) == 0) stop("Yに有効な列がありません。")

  Y <- do.call(cbind, mats)
  groups <- unlist(grps, use.names = FALSE)

  # 念のため
  if (nrow(Y) != n) stop("内部エラー: ダミー化後の行数が一致しません。")

  list(Y = Y, groups = groups)
}



guided_variance <- function(X, Y, center_X = FALSE, normalize_Y = TRUE,
                            per_feature = TRUE, drop_zero_cols = TRUE) {
  if (!is.matrix(X)) X <- as.matrix(X)

  # Y を dense ダミー行列化＋グループ名
  d <- .dummify_with_groups(Y)
  Ymat   <- d$Y
  groups <- d$groups

  # 全ゼロ列のみ除去（単一レベルfactorは残す）
  if (drop_zero_cols) {
    l2sq <- colSums(Ymat^2)
    keep <- l2sq > 0
    if (any(!keep)) {
      Ymat   <- Ymat[, keep, drop = FALSE]
      groups <- groups[keep]
    }
    if (ncol(Ymat) == 0) stop("Yの有効列がすべてゼロ列でした。")
  }

  # X 列中心化
  if (center_X) X <- scale(X, center = TRUE, scale = FALSE)

  # Y 列を L2=1 に正規化
  if (normalize_Y) {
    ynorm <- sqrt(colSums(Ymat^2))
    ynorm[ynorm == 0] <- 1
    Ymat  <- sweep(Ymat, 2, ynorm, "/")
  }

  # 計算
  XtY <- crossprod(X, Ymat)       # p × q
  num  <- sum(XtY * XtY)          # = ||X^T Y||_F^2
  den  <- sum(X * X)              # = ||X||_F^2
  frac <- if (den > 0) num / den else NA_real_

  if (!per_feature) return(frac)

  contrib_dummy <- colSums(XtY^2)
  if (!is.null(colnames(Ymat))) names(contrib_dummy) <- colnames(Ymat)

  contrib_grouped <- tapply(contrib_dummy, groups, sum)
  contrib_grouped_frac <- contrib_grouped / den  # 0–1 の比率

  list(
    fraction             = frac,                 # 全体の比率
    numerator            = num,
    denominator          = den,
    contrib              = contrib_dummy,        # ダミー列ごとの寄与
    contrib_grouped      = contrib_grouped,      # 元列ごとの寄与（生値）
    contrib_grouped_frac = contrib_grouped_frac  # 元列ごとの比率（0–1）
  )
}

# Guided PCA
results <- guided_variance(X, Y, center_X = TRUE, normalize_Y = TRUE, per_feature = TRUE)

# Save
save(results, file = outfile)