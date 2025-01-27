library("akima")
library("igraph")

args <- commandArgs(trailingOnly = TRUE)
outfile1 <- args[1]
outfile2 <- args[2]
outfile3 <- args[3]
outfile4 <- args[4]

# Load
E_cont <- unlist(read.table("plot/hpbase/cont/Landscaper/E.tsv", header=FALSE))
E_DAPT <- unlist(read.table("plot/hpbase/DAPT/Landscaper/E.tsv", header=FALSE))
E_cont_cov <- unlist(read.table("plot/hpbase/cont_cov/Landscaper/E.tsv", header=FALSE))
E_DAPT_cov <- unlist(read.table("plot/hpbase/DAPT_cov/Landscaper/E.tsv", header=FALSE))

G_sub_cont <- unlist(read.table("plot/hpbase/cont/Landscaper/SubGraph.tsv", header=FALSE))
G_sub_DAPT <- unlist(read.table("plot/hpbase/DAPT/Landscaper/SubGraph.tsv", header=FALSE))
G_sub_cont_cov <- unlist(read.table("plot/hpbase/cont_cov/Landscaper/SubGraph.tsv", header=FALSE))
G_sub_DAPT_cov <- unlist(read.table("plot/hpbase/DAPT_cov/Landscaper/SubGraph.tsv", header=FALSE))

Basin_cont <- unlist(read.table("plot/hpbase/cont/Landscaper/Basin.tsv", header=FALSE))
Basin_DAPT <- unlist(read.table("plot/hpbase/DAPT/Landscaper/Basin.tsv", header=FALSE))
Basin_cont_cov <- unlist(read.table("plot/hpbase/cont_cov/Landscaper/Basin.tsv", header=FALSE))
Basin_DAPT_cov <- unlist(read.table("plot/hpbase/DAPT_cov/Landscaper/Basin.tsv", header=FALSE))

Coordinate <- as.matrix(read.table("plot/hpbase/integrated/Landscaper/Coordinate.tsv", header=FALSE))
load("plot/hpbase/integrated/Landscaper/igraph.RData")

# Preprocess
tmp_cont <- interp(Coordinate[, 1], Coordinate[, 2], E_cont, nx = 200, ny = 200)
tmp_DAPT <- interp(Coordinate[, 1], Coordinate[, 2], E_DAPT, nx = 200, ny = 200)
tmp_cont_cov <- interp(Coordinate[, 1], Coordinate[, 2], E_cont_cov, nx = 200, ny = 200)
tmp_DAPT_cov <- interp(Coordinate[, 1], Coordinate[, 2], E_DAPT_cov, nx = 200, ny = 200)

zlim_cont <- tmp_cont$z |> as.vector() |> na.omit() |> range()
zlim_DAPT <- tmp_DAPT$z |> as.vector() |> na.omit() |> range()
zlim_cont_cov <- tmp_cont_cov$z |> as.vector() |> na.omit() |> range()
zlim_DAPT_cov <- tmp_DAPT_cov$z |> as.vector() |> na.omit() |> range()
zlim <- range(c(zlim_cont, zlim_DAPT, zlim_cont_cov, zlim_DAPT_cov))

# Plot
png(file=outfile1, width=1000, height=1000)
filled.contour(tmp_cont, color.palette = colorRampPalette(topo.colors(11, alpha = 1)), xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main="Energy"),
    axes = FALSE,
    plot.axes = {
    contour(tmp_cont$x, tmp_cont$y, tmp_cont$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.label.color = "black",
        vertex.label.cex = 0.75,
        vertex.size = 1.5,
        vertex.color = factor(G_sub_cont),
        edge.arrow.size = 0.1, add = TRUE, axes=FALSE)
    })
dev.off()

png(file=outfile2, width=1000, height=1000)
filled.contour(tmp_DAPT, color.palette = colorRampPalette(topo.colors(11, alpha = 1)), xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main="Energy"),
    axes = FALSE,
    plot.axes = {
    contour(tmp_DAPT$x, tmp_DAPT$y, tmp_DAPT$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.label.color = "black",
        vertex.label.cex = 0.75,
        vertex.size = 1.5,
        vertex.color = factor(G_sub_cont),
        edge.arrow.size = 0.1, add = TRUE, axes=FALSE)
    })
dev.off()

png(file=outfile3, width=1000, height=1000)
filled.contour(tmp_cont_cov, color.palette = colorRampPalette(topo.colors(11, alpha = 1)), xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main="Energy"),
    axes = FALSE,
    plot.axes = {
    contour(tmp_cont_cov$x, tmp_cont_cov$y, tmp_cont_cov$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.label.color = "black",
        vertex.label.cex = 0.75,
        vertex.size = 1.5,
        vertex.color = factor(G_sub_cont),
        edge.arrow.size = 0.1, add = TRUE, axes=FALSE)
    })
dev.off()

png(file=outfile4, width=1000, height=1000)
filled.contour(tmp_DAPT_cov, color.palette = colorRampPalette(topo.colors(11, alpha = 1)), xlim = c(-1, 1), ylim = c(-1, 1), zlim = zlim, asp = 1, nlevels = 24,
    key.title = title(main="Energy"),
    axes = FALSE,
    plot.axes = {
    contour(tmp_DAPT_cov$x, tmp_DAPT_cov$y, tmp_DAPT_cov$z, add = TRUE,
        levels = seq(-10, 10, by=1), col = "gray75")
    plot(g, layout=Coordinate,
        vertex.label.color = "black",
        vertex.label.cex = 0.75,
        vertex.size = 1.5,
        vertex.color = factor(G_sub_cont),
        edge.arrow.size = 0.1, add = TRUE, axes=FALSE)
    })
dev.off()
