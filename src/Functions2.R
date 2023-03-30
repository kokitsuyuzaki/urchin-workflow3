library("Seurat")
library("SeuratDisk")
library("velocyto.R")
library("SeuratWrappers")
library("ggplot2")
library("scran")
library("scater")

sample_names <- c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

group_names <- c('cont-24h', 'cont-36h', 'cont-48h', 'cont-72h', 'cont-96h', 'DAPT-24h', 'DAPT-36h', 'DAPT-48h', 'DAPT-72h', 'DAPT-96h')

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