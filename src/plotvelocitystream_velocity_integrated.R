source("src/Functions2.R")

# Parameter
db <- commandArgs(trailingOnly=TRUE)[1]
mode <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
veloList <- lapply(sample_names, function(x){
    infile1 <- paste0("output/", db, "/", x, "/velociraptor_", mode, ".RData")
    load(infile1)
    velo.obj
})

infile2 <- paste0("output/", db, "/integrated/seurat.RData")
load(infile2)

seuratList <- .stratifySeurat(seurat.integrated, group_names)

gList <- list()
length(gList) <- length(seuratList)

for(i in seq_along(group_names)){
    # Embed Velocity on UMAP coordinates
    umap.coordinates <- seuratList[[i]]@reductions$umap@cell.embeddings
    rownames(umap.coordinates) <- colnames(veloList[[i]])
    colnames(umap.coordinates) <- c("V1", "V2")
    reducedDims(veloList[[i]]) <- SimpleList(UMAP = umap.coordinates)

    # Embedding
    em <- embedVelocity(reducedDim(veloList[[i]], "UMAP"), veloList[[i]])

    # Plot object
    if(length(which(is.nan(em))) == 0){
        gList[[i]] <- plotVelocityStream(veloList[[i]], em, color.streamlines=TRUE, use.dimred="UMAP") + labs(title = group_names[[i]]) + xlim(c(-15,15)) + ylim(c(-15,15))
    }else{
        gList[[i]] <- ggplot()
    }
}
names(gList) <- group_names
gList <- gList[order(names(gList))]

# Plot
png(file=outfile, width=2000, height=1000)
(gList[[1]] | gList[[2]] | gList[[3]] | gList[[4]]) / (gList[[5]] | gList[[6]] | gList[[7]] | gList[[8]])
dev.off()
