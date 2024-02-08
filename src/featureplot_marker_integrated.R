source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
load(infile2)

marker <- marker$GENENAME_Hp
marker <- intersect(marker, rownames(seurat.integrated))

seuratList <- .stratifySeurat(seurat.integrated, group_names)

for(i in seq_along(marker)){
    filename1 <- paste0(gsub("FINISH_marker", "", outfile), marker[i], ".png")
    # replace "/" in the gene name with "_"
    if(length(grep("Hp-.*/.*.png", filename1)) != 0){
        front <- gsub("FINISH_marker", "", outfile)
        m <- regexpr("Hp-.*/.*.png", filename1)
        rear <- substr(filename1, m, m + attr(m, "match.length") - 1)
        rear <- gsub("/", "_", rear)
        filename1 <- paste0(front, rear)
    }
    if(length(grep("Sp-.*/.*.png", filename1)) != 0){
        front <- gsub("FINISH_marker", "", outfile)
        m <- regexpr("Sp-.*/.*.png", filename1)
        rear <- substr(filename1, m, m + attr(m, "match.length") - 1)
        rear <- gsub("/", "_", rear)
        filename1 <- paste0(front, rear)
    }
    # Plot
    png(file=filename1, width=600, height=600)
    print(FeaturePlot(seurat.integrated, features=marker[i], pt.size=2, label.size=6) + xlim(c(-15,15)) + ylim(c(-15,15)))
    dev.off()
    # Plot
    filename2 <- gsub(".png", "_splitby.png", filename1)
    png(file=filename2, width=2000, height=1000)
    print(.panelPlot(seuratList, group_names, marker[i]))
    dev.off()
}

# Save
file.create(outfile)
