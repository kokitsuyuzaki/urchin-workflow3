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
gene_names <- colnames(X)

# メタデータ
meta <- obj@meta.data

# ---- 予測対象（質問で挙げた列）----
targets <- c("seurat_clusters", "celltype", "germlayer", "time", "condition",
             "nFeature_SCT", "nCount_SCT", "percent.mt", "percent.rb","Phase", "dbl.dens")

# 存在チェック
targets <- targets[targets %in% colnames(meta)]
if (length(targets) == 0) stop("指定したmeta.data列が見つかりません。")

set.seed(42)

# ---- 共通指標付きCV ----
eval_one_target <- function(y_raw, X, target_name, nfolds = 5, ntree = 500) {
  is_reg <- is.numeric(y_raw) || is.integer(y_raw)
  if (!is_reg) {
    y <- factor(y_raw)
  } else {
    y <- as.numeric(y_raw)
  }

  keep <- !is.na(y) & is.finite(y)
  Xk <- X[keep, , drop = FALSE]
  yk <- y[keep]
  n  <- length(yk)

  folds <- caret::createFolds(yk, k = nfolds, list = TRUE, returnTrain = FALSE)

  # 累積器（共通）
  total_n <- 0L
  mse_sum <- 0.0           # 回帰MSE or 分類Brier（foldで合計後 / n）
  mse0_sum <- 0.0          # ベースライン（回帰=訓練平均、分類=訓練のクラス頻度）
  acc_right <- 0L          # 分類の正解数
  acc_total <- 0L

  if (!is_reg) {
    levels_all <- levels(yk)  # 全データでのクラス集合
  }

  for (i in seq_along(folds)) {
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(n), test_idx)

    Xtr <- Xk[train_idx, , drop = FALSE]
    Xte <- Xk[test_idx, , drop = FALSE]
    ytr <- yk[train_idx]
    yte <- yk[test_idx]
    nte <- length(yte)
    total_n <- total_n + nte

    if (is_reg) {
      # --- 回帰 ---
      fit <- ranger(x = as.data.frame(Xtr), y = as.numeric(ytr),
                    num.trees = ntree, mtry = max(1, floor(sqrt(ncol(Xtr)))),
                    min.node.size = 5, seed = 42,
                    classification = FALSE, keep.inbag = FALSE)
      pred <- as.numeric(predict(fit, data = as.data.frame(Xte))$predictions)

      # MSE（共通指標）
      mse_sum  <- mse_sum  + sum((as.numeric(yte) - pred)^2)
      # ベースライン：訓練平均
      mu_tr    <- mean(as.numeric(ytr))
      mse0_sum <- mse0_sum + sum((as.numeric(yte) - mu_tr)^2)

    } else {
      # --- 分類 ---
      # 確率を出す（Brier用）
      fit <- ranger(x = as.data.frame(Xtr), y = factor(ytr, levels = levels_all),
                    num.trees = ntree, mtry = max(1, floor(sqrt(ncol(Xtr)))),
                    min.node.size = 5, seed = 42,
                    classification = TRUE, probability = TRUE)
      pr <- predict(fit, data = as.data.frame(Xte))$predictions
      # 列（クラス）を全クラスに整列
      probs <- matrix(0, nrow = nte, ncol = length(levels_all))
      colnames(probs) <- levels_all
      common <- intersect(colnames(pr), levels_all)
      probs[, match(common, levels_all)] <- pr[, common, drop = FALSE]

      # 予測ラベルとAccuracy
      yhat <- factor(colnames(probs)[max.col(probs, ties.method = "first")], levels = levels_all)
      acc_right <- acc_right + sum(yhat == factor(yte, levels = levels_all))
      acc_total <- acc_total + nte

      # One-hot 真値
      Yte <- diag(length(levels_all))[as.integer(factor(yte, levels = levels_all)), , drop = FALSE]
      # Brier（共通MSE）
      mse_sum <- mse_sum + sum(rowSums((probs - Yte)^2))

      # ベースライン：訓練データのクラス頻度（climatology）
      p_bar <- prop.table(table(factor(ytr, levels = levels_all)))
      P0 <- matrix(rep(as.numeric(p_bar), each = nte), nrow = nte, ncol = length(levels_all))
      mse0_sum <- mse0_sum + sum(rowSums((P0 - Yte)^2))
    }
  }

  # 集計
  common_mse   <- mse_sum / total_n
  common_mse0  <- mse0_sum / total_n
  common_skill <- 1 - common_mse / common_mse0  # 回帰=R2, 分類=Brier Skill Score

  if (is_reg) {
    # 参考: R2, RMSE（全fold結合の定義に合わせる）
    R2   <- common_skill
    RMSE <- sqrt(common_mse)
    tibble(
      target = target_name, task = "regression", n = total_n,
      common_MSE = common_mse, common_baseline = common_mse0, common_skill = R2,
      R2 = R2, RMSE = RMSE
    )
  } else {
    Accuracy <- acc_right / acc_total
    tibble(
      target = target_name, task = "classification", n = total_n, n_class = length(levels_all),
      common_MSE = common_mse, common_baseline = common_mse0, common_skill = common_skill,
      Accuracy = Accuracy
    )
  }
}

# ---- 実行（木本数は重いなら下げる） ----
results <- map_dfr(targets, ~ eval_one_target(meta[[.x]], X, .x, nfolds = 5, ntree = 100))

# 共通指標で並べ替え
results <- results %>% arrange(desc(common_skill))

# 保存
save(results, file = outfile)
