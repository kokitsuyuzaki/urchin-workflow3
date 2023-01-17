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

markers <- c("Hp-Hnf6", "Hp-Chordin", "Hp-Hox7", "Hp-FoxJ1", "Hp-BMP2_4", "Hp-Bra", "Hp-Delta", "Hp-SoxC", "Hp-Awh", "Hp-Gcm", "Hp-Blimp1", "Hp-Endo16", "Hp-Ephrin", "Hp-IrxA", "Hp-Hox11_13b", "Hp-FoxQ2", "Hp-MSP130", "Hp-Erg", "Hp-Ese", "Hp-FoxY", "Hp-Echn38", "Hp-Tph", "Hp-Pnlip-5", "Hp-Alx1", "Hp-Sm50")

sample_colors <- brewer.pal(8, "Dark2")

sample_names <- c('cont-24h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

group_names <- c('cont-24h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

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
    list(seurat.each1, seurat.each2, seurat.each3, seurat.each4,
        seurat.each5, seurat.each6, seurat.each7, seurat.each8)
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
    (gList[[1]] | gList[[2]] | gList[[3]] | gList[[4]]) / (gList[[5]] | gList[[6]] | gList[[7]] | gList[[8]])
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
    cont_data <- data[data[, "sample"] %in% c('cont-24h', 'cont-48h', 'cont-72h', 'cont-96h'), ]
    DAPT_data <- data[data[, "sample"] %in% c('DAPT-24h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h'), ]
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