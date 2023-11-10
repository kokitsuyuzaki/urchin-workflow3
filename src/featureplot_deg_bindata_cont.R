source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
numSheets <- length(excel_sheets(infile2))

# Preprocessing
dir.create(gsub("FINISH", "", outfile))

for(i in seq_len(numSheets)){
    # Loading
    marker <- read_excel(infile2, i)
    if(nrow(marker) != 0){
        # Preprocessing
        marker <- intersect(unlist(marker[,1]), rownames(seurat.cont))
        out.dir <- paste0(gsub("FINISH", "", outfile), "variable", i)
        dir.create(out.dir)
        # Plot
        for(j in seq_len(min(length(marker), 50))){
            filename <- paste0(out.dir, "/", marker[j], ".png")
            # replace "/" in the gene name with "_"
            if(length(grep("Hp-.*/.*.png", filename)) != 0){
                front <- gsub("FINISH", "", outfile)
                m <- regexpr("Hp-.*/.*.png", filename)
                rear <- substr(filename, m, m + attr(m, "match.length") - 1)
                rear <- gsub("/", "_", rear)
                filename <- paste0(front, rear)
            }
            # Plot
            png(file=filename, width=600, height=600)
            print(FeaturePlot(seurat.cont, features=marker[j], pt.size=2, label.size=6))
            dev.off()
        }
    }
}

# Save
file.create(outfile)