library("GSEABase")
library("GO.db")
library("Seurat")
library("sctransform")
library("readxl")
library("writexl")
library("scDblFinder")
library("scran")
library("scater")
library("SeuratWrappers")
library("monocle3")
library("velociraptor")
library("scTGIF")
library("RColorBrewer")
library("reticulate")
library("abind")
library("purrr")
library("tidyr")
library("stringr")
library("ggplot2")
library("viridis")
library("reshape2")
library("qvalue")

# Label Stratification
.labelStratify <- function(seurat.integrated, target.pattern){
    target <- grep(target.pattern, seurat.integrated@meta.data$sample)
    seurat.integrated[, target]
}

germlayer_colors <- c(brewer.pal(8, "Dark2")[c(2,3,1)], rgb(0,0,0, 0.5))
names(germlayer_colors) <- c("Endoderm", "Mesoderm", "Ectoderm", "NA")

# rm Hp-BMP2_4, Hp-Hox11_13b, Hp-FoxQ2, Hp-MSP130
markers <- c("Hp-Hnf6", "Hp-Chordin", "Hp-Hox7", "Hp-FoxJ1", "Hp-Bra", "Hp-Delta", "Hp-SoxC", "Hp-Awh", "Hp-Gcm", "Hp-Blimp1", "Hp-Endo16", "Hp-Ephrin", "Hp-IrxA", "Hp-Erg", "Hp-Ese", "Hp-FoxY", "Hp-Echn38", "Hp-Tph", "Hp-Pnlip-5", "Hp-Alx1", "Hp-Sm50")

fig2_markers <- c("Hp-SoxC", "Hp-Hnf6", "Hp-FoxJ1", "Hp-Chordin", "Hp-Ephrin", "Hp-Hox7", "Hp-Gcm", "Hp-Echn38", "Hp-FoxY", "Hp-Erg", "Hp-Ese", "Hp-Sm50", "Hp-Alx1", "Hp-Endo16", "Hp-Blimp1", "Hp-Bra", "Hp-IrxA", "Hp-Pnlip-5")

fig4_markers <- c("Hp-Gad", "Hp-Hdc-3", "Hp-Th", "Hp-Chat", "Hp-Ddc", "Hp-Tph", "Hp-Syt1-1-like", "Hp-Otp", "Hp-Ngn", "Hp-Ac/Sc", "Hp-Awh", "Hp-Delta", "Hp-SoxC", "Hp-Smad-ip", "Hp-Rx", "Hp-Hbn", "Hp-Nkx2.1", "Hp-Six3", "Hp-FoxQ2-1-like", "Hp-FoxQ2-1", "Hp-HesC", "Hp-Soxb1")

fig2_celltypes <- c("uncharacterized", "Pancreas", "Anus", "Endoderm", "Stomach_Intestine", "Intestine", "Stomach", "Skeleton", "Non_skeleton_mesoderm", "Blastocoelar_cell", "Germ_line_future", "Pigment", "Aboral_ectoderm", "Oral_ectoderm", "Ciliary_band", "Neurons")

sample_colors <- c(brewer.pal(9, "Set1"), "black")

sample_names <- c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

group_names <- c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

conditions <- c("cont", "DAPT")

times <- c("24h", "36h", "48h", "72h", "96h")

seurat_clusters <- c(
    "3", "29", "12", "1", "17", "41", "33", "23", "38", "40",
    "25", "32", "36", "35", "39", "21", "30", "13", "43", "27",
    "14", "4", "42", "2", "19", "31", "9", "37", "15", "22",
    "8", "34", "11", "0", "16", "26", "28", "6", "5", "24",
    "20", "18", "7", "10")

blacklist_gene <- "Sp-Hrh2_3"

genes_ridgeplot_hpbase <- c("Hp-Pcna", "Hp-Srrm2-like", "Hp-Tpx2L1",
    "Hp-Lbr", "Hp-Map215prh-like", "Hp-Ndc80L")

genes_ridgeplot_echinobase <- c("Sp-Pcna", "Sp-Srrm2", "Sp-Tpx2L1",
    "Sp-Map215prh", "Sp-Ndc80L")

.firstup <- function(x){
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

.tableToGMT <- function(gotable){
    # GeneSet's List
    gsc <- lapply(unique(gotable$TERM), function(g){
        target <- which(gotable$TERM == g)
        geneIds <- gotable$GENENAME[target]
        GeneSet(setName=g, geneIds)
    })
    GeneSetCollection(gsc)
}

.stratifySeurat <- function(seurat.integrated, group_names){
    seurat.each1 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[1])]
    seurat.each2 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[2])]
    seurat.each3 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[3])]
    seurat.each4 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[4])]
    seurat.each5 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[5])]
    seurat.each6 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[6])]
    seurat.each7 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[7])]
    seurat.each8 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[8])]
    seurat.each9 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[9])]
    seurat.each10 <- seurat.integrated[,
        which(seurat.integrated@meta.data$sample == group_names[10])]
    list(seurat.each1, seurat.each2, seurat.each3, seurat.each4,
        seurat.each5, seurat.each6, seurat.each7, seurat.each8, seurat.each9, seurat.each10)
}

.panelPlot <- function(seuratList, group_names, features){
    # Stratify
    gList <- list()
    length(gList) <- length(group_names)
    for(i in seq_along(group_names)){
        gList[[i]] <- FeaturePlot(seuratList[[i]], features=features,
            reduction = "umap", pt.size=2, label.size=6) + labs(title=group_names[i]) + xlim(c(-15,15)) + ylim(c(-15,15))

    }
    names(gList) <- group_names
    gList <- gList[order(names(gList))]
    # Patch work
    (gList[[1]] | gList[[2]] | gList[[3]] | gList[[4]] | gList[[5]]) / (gList[[6]] | gList[[7]] | gList[[8]] | gList[[9]] | gList[[10]])
}

.BarPlot <- function(seurat.integrated){
    data <- seurat.integrated@meta.data[,
        c("seurat_clusters", "sample")]
    data <- .countCategory(data)
    data <- .sortByDiff(data)
    g <- ggplot(data, aes(x=seurat_clusters, y=no_cells, fill=seurat_clusters))
    g <- g + geom_bar(stat='identity')
    g <- g + xlab('Cell type')
    g <- g + ylab('# cells')
    g <- g + facet_wrap(~sample, ncol=1)
    g <- g + theme(strip.text.x = element_text(size = 25))
    g
}

.countCategory <- function(data){
    x <- as.character(unique(data[, "seurat_clusters"]))
    y <- as.character(unique(data[, "sample"]))
    xy <- as.matrix(expand.grid(x, y))
    out <- do.call("rbind", apply(xy, 1, function(z){
        targetx <- which(data[, "seurat_clusters"] == z[1])
        targety <- which(data[, "sample"] == z[2])
        no_cells <- length(intersect(targetx, targety))
        data.frame(z[1], z[2], no_cells)
    }))
    colnames(out) <- c("seurat_clusters", "sample", "no_cells")
    out
}

.sortByDiff <- function(data){
    cont_data <- data[data[, "sample"] %in% c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h'), ]
    DAPT_data <- data[data[, "sample"] %in% c('DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h'), ]
    sum_cont_data <- .sumByCelltype(cont_data)
    sum_DAPT_data <- .sumByCelltype(DAPT_data)
    diff <- sum_DAPT_data[,2] - sum_cont_data[,2]
    sum_DAPT_data[,1][order(diff)]
    data[, "seurat_clusters"] <- factor(data[, "seurat_clusters"], levels=sum_DAPT_data[,1][order(diff)])
    data
}

.sumByCelltype <- function(data){
    target <- unique(data[,1])
    no_cells <- do.call("rbind", lapply(target, function(x){
        sum(data[which(data[,1] == x), "no_cells"])
    }))
    data.frame(seurat_clusters=target, no_cells=no_cells)
}
