source("src/Functions.R")

# Loading
load("output/hpbase/integrated/seurat.RData")

seuratList <- .stratifySeurat(seurat.integrated, group_names)

genes <- c("Hp-Opn4L", "Hp-Opn4L-1", "Hp-Opn5L", "Hp-Opsin1", "Hp-Opsin4", "Hp-Opsin1-like", "Hp-Opsin2", "Hp-Opsin5-like", "Hp-Opsin3.2-like", "Hp-Nos-2", "Hp-Nos1")

for(i in seq_along(genes)){
    filename = paste0(genes[i], ".png")
    png(file=filename, width=2000, height=1000)
    plot(.panelPlot(seuratList, group_names, genes[i]))
    dev.off()
}


