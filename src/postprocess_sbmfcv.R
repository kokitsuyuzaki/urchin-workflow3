outfile <- commandArgs(trailingOnly=TRUE)[1]

# Preprpcessing
nmf_rank <- ncol(read.table('output/hpbase/integrated/sbmfcv/BIN_DATA.tsv', header=FALSE, nrow=1))

out <- lapply(seq(50),
    function(i){
        load(paste0("output/hpbase/integrated/sbmfcv/nmf/", nmf_rank, "/", i, ".RData"))
        out$U})
out_mat <- do.call(cbind, out)

# k-meansクラスタリング
set.seed(1234)  # 再現性確保
km <- kmeans(t(out_mat), centers = nmf_rank, nstart = 50)  # 列をクラスタリング
# km$cluster: 各列のクラスタ番号（1..rank_nmf）
# km$centers: 各クラスタの平均基底ベクトル（列方向の平均）

# クラスタ平均を「平均的な基底」として再構成
U_consensus <- t(km$centers)

# 2値化
U_consensus[which(U_consensus < 0.5)] <- -1
U_consensus[which(U_consensus >= 0.5)] <- 1

# Save
write.table(U_consensus, outfile, row.names=FALSE, col.names=FALSE, quote=FALSE)