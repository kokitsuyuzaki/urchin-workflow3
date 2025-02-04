# Package Loading
library("dynwrap")
library("dyno")
library("dynmethods")
library("babelwhale")
library("Seurat")
library("tidyverse")

# Singularityの設定
test_singularity_installation(detailed = TRUE)
options(dynwrap_backend = 'container')

# Cacheの設定
config <- create_singularity_config(cache_dir = "./cache_dir")
set_default_config(config)

# クラスタ中心に最近傍データの取得
.nearest_to_centroid <- function(mat, cluster_col){
  # クラスタごとの中心を計算
  id <- rownames(mat)
  data <- data.frame(id = id, cluster=cluster_col, as.data.frame(mat))
  centroids <- data %>%
    group_by(cluster) %>%
    summarise(across(where(is.numeric), mean), .groups = "drop")

  # クラスタごとに最近傍データを取得
  nearest_ids <- data %>%
    group_by(cluster) %>%
    mutate(
        distance = rowSums((across(where(is.numeric)) - centroids[match(cluster, centroids$cluster), -1])^2)^0.5
    ) %>%
    slice_min(distance) %>%
    select(id)
  nearest_ids$id
}