source("src/Functions.R")

# Parameter
db <- commandArgs(trailingOnly=TRUE)[1]
infile1 <- commandArgs(trailingOnly=TRUE)[2]
infile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile <- commandArgs(trailingOnly=TRUE)[4]

# Loading
load(infile1)
load(infile2)

if(db == "hpbase"){
    marker <- marker$GENENAME_Hp
}else{
    marker <- marker$GENENAME_Sp
}
marker <- intersect(marker, rownames(seurat.obj))

for(i in seq_along(marker)){
    filename <- paste0(gsub("FINISH_marker", "", outfile), marker[i], ".png")
    # replace "/" in the gene name with "_"
    if(length(grep("Hp-.*/.*.png", filename)) != 0){
        front <- gsub("FINISH_marker", "", outfile)
        m <- regexpr("Hp-.*/.*.png", filename)
        rear <- substr(filename, m, m + attr(m, "match.length") - 1)
        rear <- gsub("/", "_", rear)
        filename <- paste0(front, rear)
    }
    if(length(grep("Sp-.*/.*.png", filename)) != 0){
        front <- gsub("FINISH_marker", "", outfile)
        m <- regexpr("Sp-.*/.*.png", filename)
        rear <- substr(filename, m, m + attr(m, "match.length") - 1)
        rear <- gsub("/", "_", rear)
        filename <- paste0(front, rear)
    }
    # Plot
    png(file=filename, width=600, height=600)
    print(FeaturePlot(seurat.obj, features=marker[i], pt.size=2, label.size=6))
    dev.off()
}

# Save
file.create(outfile)