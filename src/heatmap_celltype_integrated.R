source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]
outfile4 <- commandArgs(trailingOnly=TRUE)[5]
outfile5 <- commandArgs(trailingOnly=TRUE)[6]

# Loading
load(infile)

# Preprocessing
data <- seurat.integrated[fig3_markers, ]

cont24h <- as.matrix(data[, which(data@meta.data$sample == "cont-24h")]@assays$RNA@data)
cont36h <- as.matrix(data[, which(data@meta.data$sample == "cont-36h")]@assays$RNA@data)
cont48h <- as.matrix(data[, which(data@meta.data$sample == "cont-48h")]@assays$RNA@data)
cont72h <- as.matrix(data[, which(data@meta.data$sample == "cont-72h")]@assays$RNA@data)
cont96h <- as.matrix(data[, which(data@meta.data$sample == "cont-96h")]@assays$RNA@data)
cont_ <- as.matrix(data[, grep("cont-", data@meta.data$sample)]@assays$RNA@data)

DAPT24h <- as.matrix(data[, which(data@meta.data$sample == "DAPT-24h")]@assays$RNA@data)
DAPT36h <- as.matrix(data[, which(data@meta.data$sample == "DAPT-36h")]@assays$RNA@data)
DAPT48h <- as.matrix(data[, which(data@meta.data$sample == "DAPT-48h")]@assays$RNA@data)
DAPT72h <- as.matrix(data[, which(data@meta.data$sample == "DAPT-72h")]@assays$RNA@data)
DAPT96h <- as.matrix(data[, which(data@meta.data$sample == "DAPT-96h")]@assays$RNA@data)
DAPT_ <- as.matrix(data[, grep("DAPT-", data@meta.data$sample)]@assays$RNA@data)

cont <- cbind(
	rowMeans(cont24h),
	rowMeans(cont36h),
	rowMeans(cont48h),
	rowMeans(cont72h),
	rowMeans(cont96h),
	rowMeans(cont_))
cont <- log10(cont)
colnames(cont) <- c(times, "overall")

DAPT <- cbind(
	rowMeans(DAPT24h),
	rowMeans(DAPT36h),
	rowMeans(DAPT48h),
	rowMeans(DAPT72h),
	rowMeans(DAPT96h),
	rowMeans(DAPT_))
DAPT <- log10(DAPT)
colnames(DAPT) <- c(times, "overall")

Ratio <- DAPT - cont

T_Pval24h <- sapply(seq_along(fig3_markers), function(x){
	t.test(cont24h[x, ], DAPT24h[x, ])$p.value
})
T_Pval36h <- sapply(seq_along(fig3_markers), function(x){
	t.test(cont36h[x, ], DAPT36h[x, ])$p.value
})
T_Pval48h <- sapply(seq_along(fig3_markers), function(x){
	t.test(cont48h[x, ], DAPT48h[x, ])$p.value
})
T_Pval72h <- sapply(seq_along(fig3_markers), function(x){
	t.test(cont72h[x, ], DAPT72h[x, ])$p.value
})
T_Pval96h <- sapply(seq_along(fig3_markers), function(x){
	t.test(cont96h[x, ], DAPT96h[x, ])$p.value
})
T_Pval_ <- sapply(seq_along(fig3_markers), function(x){
	t.test(cont_[x, ], DAPT_[x, ])$p.value
})

T_Pvalue <- cbind(
	T_Pval24h,
	T_Pval36h,
	T_Pval48h,
	T_Pval72h,
	T_Pval96h,
	T_Pval_)
rownames(T_Pvalue) <- fig3_markers
colnames(T_Pvalue) <- c(times, "overall")
# T_Pvalue <- -log10(T_Pvalue)

W_Pval24h <- sapply(seq_along(fig3_markers), function(x){
	wilcox.test(cont24h[x, ], DAPT24h[x, ])$p.value
})
W_Pval36h <- sapply(seq_along(fig3_markers), function(x){
	wilcox.test(cont36h[x, ], DAPT36h[x, ])$p.value
})
W_Pval48h <- sapply(seq_along(fig3_markers), function(x){
	wilcox.test(cont48h[x, ], DAPT48h[x, ])$p.value
})
W_Pval72h <- sapply(seq_along(fig3_markers), function(x){
	wilcox.test(cont72h[x, ], DAPT72h[x, ])$p.value
})
W_Pval96h <- sapply(seq_along(fig3_markers), function(x){
	wilcox.test(cont96h[x, ], DAPT96h[x, ])$p.value
})
W_Pval_ <- sapply(seq_along(fig3_markers), function(x){
	wilcox.test(cont_[x, ], DAPT_[x, ])$p.value
})

W_Pvalue <- cbind(
	W_Pval24h,
	W_Pval36h,
	W_Pval48h,
	W_Pval72h,
	W_Pval96h,
	W_Pval_)
rownames(W_Pvalue) <- fig3_markers
colnames(W_Pvalue) <- c(times, "overall")
# W_Pvalue <- -log10(W_Pvalue)

# Plot
g_cont <- melt(cont)
colnames(g_cont) <- c("Symbol", "Time", "Value")
g_cont <- ggplot(g_cont, aes(x=Time, y=Symbol, fill=Value))
g_cont <- g_cont + geom_tile()
g_cont <- g_cont + scale_fill_gradientn(colors = viridis(256), name = "Log10(Avg. Exp.)")
g_cont <- g_cont + labs(x=NULL, y=NULL)
ggsave(file=outfile1, g_cont, dpi=200, width=6, height=7)

g_DAPT <- melt(DAPT)
colnames(g_DAPT) <- c("Symbol", "Time", "Value")
g_DAPT <- ggplot(g_DAPT, aes(x=Time, y=Symbol, fill=Value))
g_DAPT <- g_DAPT + geom_tile()
g_DAPT <- g_DAPT + scale_fill_gradientn(colors = viridis(256), name = "Log10(Avg. Exp.)")
g_DAPT <- g_DAPT + labs(x=NULL, y=NULL)
ggsave(file=outfile2, g_DAPT, dpi=200, width=6, height=7)

g_Ratio <- melt(Ratio)
colnames(g_Ratio) <- c("Symbol", "Time", "Value")
g_Ratio <- ggplot(g_Ratio, aes(x=Time, y=Symbol, fill=Value))
g_Ratio <- g_Ratio + geom_tile()
g_Ratio <- g_Ratio + scale_fill_gradientn(colors = viridis(256), name = "Log10(Ratio)")
g_Ratio <- g_Ratio + labs(x=NULL, y=NULL)
ggsave(file=outfile3, g_Ratio, dpi=200, width=6, height=7)

g_T_Pvalue <- melt(T_Pvalue)
colnames(g_T_Pvalue) <- c("Symbol", "Time", "Value")
g_T_Pvalue <- ggplot(g_T_Pvalue, aes(x=Time, y=Symbol, fill=Value))
g_T_Pvalue <- g_T_Pvalue + geom_tile()
g_T_Pvalue <- g_T_Pvalue + scale_fill_gradientn(colors = viridis(256), name = "Pvalue")
g_T_Pvalue <- g_T_Pvalue + labs(x=NULL, y=NULL)
ggsave(file=outfile4, g_T_Pvalue, dpi=200, width=6, height=7)

g_W_Pvalue <- melt(W_Pvalue)
colnames(g_W_Pvalue) <- c("Symbol", "Time", "Value")
g_W_Pvalue <- ggplot(g_W_Pvalue, aes(x=Time, y=Symbol, fill=Value))
g_W_Pvalue <- g_W_Pvalue + geom_tile()
g_W_Pvalue <- g_W_Pvalue + scale_fill_gradientn(colors = viridis(256), name = "Pvalue")
g_W_Pvalue <- g_W_Pvalue + labs(x=NULL, y=NULL)
ggsave(file=outfile5, g_W_Pvalue, dpi=200, width=6, height=7)
