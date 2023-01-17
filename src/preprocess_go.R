source("src/Functions.R")

# Parameter
infile <- commandArgs(trailingOnly=TRUE)[1]
outfile1 <- commandArgs(trailingOnly=TRUE)[2]
outfile2 <- commandArgs(trailingOnly=TRUE)[3]
outfile3 <- commandArgs(trailingOnly=TRUE)[4]
outfile4 <- commandArgs(trailingOnly=TRUE)[5]
outfile5 <- commandArgs(trailingOnly=TRUE)[6]
outfile6 <- commandArgs(trailingOnly=TRUE)[7]

# Loading
go <- read.delim(infile, header=FALSE)

# GENENAME - GOID
genename <- unlist(apply(go, 1, function(x){
    out <- strsplit(as.character(x[4]), ",")[[1]]
    rep(x[3], length=length(out))
}))
goid <- unlist(apply(go, 1, function(x){
    strsplit(as.character(x[4]), ",")[[1]]
}))
gotable1 <- cbind(genename, goid)
colnames(gotable1) <- c("GENENAME", "GOID")

# GOID - GOTERM
gotable2 <- select(GO.db, keys=unique(gotable1[,2]), keytype="GOID",
    columns=c("GOID", "TERM", "ONTOLOGY"))

# Merge
gotable <- unique(merge(gotable1, gotable2, by="GOID"))

#############################
# Preprocessing for Hpbase
#############################
gotable_bp_hpbase <- gotable[which(gotable$ONTOLOGY == "BP"), ]
gotable_mf_hpbase <- gotable[which(gotable$ONTOLOGY == "MF"), ]
gotable_cc_hpbase <- gotable[which(gotable$ONTOLOGY == "CC"), ]

gotable_bp_hpbase <- gotable_bp_hpbase[grep("^LOC", gotable_bp_hpbase$GENENAME, invert=TRUE), ]
gotable_mf_hpbase <- gotable_mf_hpbase[grep("^LOC", gotable_mf_hpbase$GENENAME, invert=TRUE), ]
gotable_cc_hpbase <- gotable_cc_hpbase[grep("^LOC", gotable_cc_hpbase$GENENAME, invert=TRUE), ]

gotable_bp_hpbase$GENENAME <- paste0("Hp-", .firstup(tolower(gotable_bp_hpbase$GENENAME)))
gotable_mf_hpbase$GENENAME <- paste0("Hp-", .firstup(tolower(gotable_mf_hpbase$GENENAME)))
gotable_cc_hpbase$GENENAME <- paste0("Hp-", .firstup(tolower(gotable_cc_hpbase$GENENAME)))

gmt_go_bp_hpbase <- .tableToGMT(gotable_bp_hpbase)
gmt_go_mf_hpbase <- .tableToGMT(gotable_mf_hpbase)
gmt_go_cc_hpbase <- .tableToGMT(gotable_cc_hpbase)

#############################
# Preprocessing for Echinobase
#############################
gotable_bp_echinobase <- gotable[which(gotable$ONTOLOGY == "BP"), ]
gotable_mf_echinobase <- gotable[which(gotable$ONTOLOGY == "MF"), ]
gotable_cc_echinobase <- gotable[which(gotable$ONTOLOGY == "CC"), ]

gotable_bp_echinobase$GENENAME <- paste0("Sp-", .firstup(tolower(gotable_bp_echinobase$GENENAME)))
gotable_mf_echinobase$GENENAME <- paste0("Sp-", .firstup(tolower(gotable_mf_echinobase$GENENAME)))
gotable_cc_echinobase$GENENAME <- paste0("Sp-", .firstup(tolower(gotable_cc_echinobase$GENENAME)))

gmt_go_bp_echinobase <- .tableToGMT(gotable_bp_echinobase)
gmt_go_mf_echinobase <- .tableToGMT(gotable_mf_echinobase)
gmt_go_cc_echinobase <- .tableToGMT(gotable_cc_echinobase)

# Save
save(gmt_go_bp_hpbase, file=outfile1)
save(gmt_go_mf_hpbase, file=outfile2)
save(gmt_go_cc_hpbase, file=outfile3)
save(gmt_go_bp_echinobase, file=outfile4)
save(gmt_go_mf_echinobase, file=outfile5)
save(gmt_go_cc_echinobase, file=outfile6)
