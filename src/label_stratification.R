source("src/Functions.R")

# Loading
load('output/hpbase/integrated/seurat_annotated.RData')
meta.data <- seurat.integrated@meta.data

# Sample (cont)
load("output/hpbase/cont-24h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "cont-24h"), ]
save(seurat.obj, file="output/hpbase/cont-24h/seurat_annotated.RData")

load("output/hpbase/cont-36h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "cont-36h"), ]
save(seurat.obj, file="output/hpbase/cont-36h/seurat_annotated.RData")

load("output/hpbase/cont-48h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "cont-48h"), ]
save(seurat.obj, file="output/hpbase/cont-48h/seurat_annotated.RData")

load("output/hpbase/cont-72h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "cont-72h"), ]
save(seurat.obj, file="output/hpbase/cont-72h/seurat_annotated.RData")

load("output/hpbase/cont-96h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "cont-96h"), ]
save(seurat.obj, file="output/hpbase/cont-96h/seurat_annotated.RData")

# Sample (DAPT)
load("output/hpbase/DAPT-24h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "DAPT-24h"), ]
save(seurat.obj, file="output/hpbase/DAPT-24h/seurat_annotated.RData")

load("output/hpbase/DAPT-36h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "DAPT-36h"), ]
save(seurat.obj, file="output/hpbase/DAPT-36h/seurat_annotated.RData")

load("output/hpbase/DAPT-48h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "DAPT-48h"), ]
save(seurat.obj, file="output/hpbase/DAPT-48h/seurat_annotated.RData")

load("output/hpbase/DAPT-72h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "DAPT-72h"), ]
save(seurat.obj, file="output/hpbase/DAPT-72h/seurat_annotated.RData")

load("output/hpbase/DAPT-96h/seurat.RData")
seurat.obj@meta.data <- meta.data[which(meta.data$sample == "DAPT-96h"), ]
save(seurat.obj, file="output/hpbase/DAPT-96h/seurat_annotated.RData")

# Condition
load("output/hpbase/cont/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("cont", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/cont/seurat_annotated.RData")

load("output/hpbase/DAPT/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("DAPT", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/DAPT/seurat_annotated.RData")

# Time
load("output/hpbase/24h/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("24h", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/24h/seurat_annotated.RData")

load("output/hpbase/36h/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("36h", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/36h/seurat_annotated.RData")

load("output/hpbase/48h/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("48h", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/48h/seurat_annotated.RData")

load("output/hpbase/72h/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("72h", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/72h/seurat_annotated.RData")

load("output/hpbase/96h/seurat.RData")
seurat.integrated@meta.data <- meta.data[grep("96h", meta.data$sample), ]
save(seurat.integrated, file="output/hpbase/96h/seurat_annotated.RData")






# Loading for Stratification
load('output/hpbase/integrated/seurat_annotated.RData')
seurat.ref <- seurat.integrated

# Sample (cont, Stratified)
seurat.obj <- .labelStratify(seurat.ref, "cont-24h")
save(seurat.obj, file="output/hpbase/cont-24h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "cont-36h")
save(seurat.obj, file="output/hpbase/cont-36h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "cont-48h")
save(seurat.obj, file="output/hpbase/cont-48h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "cont-72h")
save(seurat.obj, file="output/hpbase/cont-72h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "cont-96h")
save(seurat.obj, file="output/hpbase/cont-96h_stratified/seurat_annotated.RData")

# Sample (DAPT, Stratified)
seurat.obj <- .labelStratify(seurat.ref, "DAPT-24h")
save(seurat.obj, file="output/hpbase/DAPT-24h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "DAPT-36h")
save(seurat.obj, file="output/hpbase/DAPT-36h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "DAPT-48h")
save(seurat.obj, file="output/hpbase/DAPT-48h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "DAPT-72h")
save(seurat.obj, file="output/hpbase/DAPT-72h_stratified/seurat_annotated.RData")

seurat.obj <- .labelStratify(seurat.ref, "DAPT-96h")
save(seurat.obj, file="output/hpbase/DAPT-96h_stratified/seurat_annotated.RData")

# Condition (Stratified)
seurat.integrated <- .labelStratify(seurat.ref, "cont-")
save(seurat.integrated, file="output/hpbase/cont_stratified/seurat_annotated.RData")

seurat.integrated <- .labelStratify(seurat.ref, "DAPT-")
save(seurat.integrated, file="output/hpbase/DAPT_stratified/seurat_annotated.RData")

# Time (Stratified)
seurat.integrated <- .labelStratify(seurat.ref, "-24h")
save(seurat.integrated, file="output/hpbase/24h_stratified/seurat_annotated.RData")

seurat.integrated <- .labelStratify(seurat.ref, "-36h")
save(seurat.integrated, file="output/hpbase/36h_stratified/seurat_annotated.RData")

seurat.integrated <- .labelStratify(seurat.ref, "-48h")
save(seurat.integrated, file="output/hpbase/48h_stratified/seurat_annotated.RData")

seurat.integrated <- .labelStratify(seurat.ref, "-72h")
save(seurat.integrated, file="output/hpbase/72h_stratified/seurat_annotated.RData")

seurat.integrated <- .labelStratify(seurat.ref, "-96h")
save(seurat.integrated, file="output/hpbase/96h_stratified/seurat_annotated.RData")
