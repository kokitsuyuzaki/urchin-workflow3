source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile <- commandArgs(trailingOnly=TRUE)[3]

# Loading
load(infile1)
bindata <- read.table(infile2, header=FALSE)
colnames(bindata) <- paste0("Var", seq(ncol(bindata)))

# Marker detection
for(i in seq_len(ncol(bindata))){
	cmd <- paste0("seurat.integrated$Var", i, "<- bindata[,", i, "]")
	eval(parse(text=cmd))
}
outList <- lapply(seq(ncol(bindata)), function(i){
     tmp <- FindMarkers(seurat.integrated, ident.1="1", group.by=paste0("Var", i))
     genenames <- rownames(tmp)
     cbind(genenames, tmp)
})
names(outList) <- colnames(bindata)

# Save
write_xlsx(outList, outfile)
