source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
h_cont <- read.table("plot/hpbase/cont/Landscaper/h.tsv", header=FALSE)
h_dapt <- read.table("plot/hpbase/dapt/Landscaper/h.tsv", header=FALSE)

# Name
rownames_h <- paste0("Var", seq(nrow(h_cont)))

# Re-ordering
h_diff <- h_dapt - h_cont
order_h_diff <- order(unlist(h_diff), decreasing=TRUE)
h_cont <- h_cont[order_h_diff, ]
h_dapt <- h_dapt[order_h_diff, ]
h_diff <- h_diff[order_h_diff, ]
rownames_h <- rownames_h[order_h_diff]

# Data for ggplot2
gdata <- data.frame(value=c(h_cont, h_dapt),
        group=c(rep("cont", length=length(h_cont)),
            rep("dapt", length=length(h_dapt))),
        variable=factor(rownames_h, level=rownames_h))

# Plot
g <- ggplot() + geom_line(aes(x=variable, y=value, colour=group, group=group), linewidth=1, data = gdata, stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(outfile, plot=g, width=7, height=4.5)
