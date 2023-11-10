source("src/Functions.R")

# Parameter
outfile <- commandArgs(trailingOnly=TRUE)[1]

# Loading
J_cont <- read.table("plot/hpbase/cont/Landscaper/J.tsv", header=FALSE)
J_dapt <- read.table("plot/hpbase/dapt/Landscaper/J.tsv", header=FALSE)

# Name
upt_J <- upper.tri(J_cont)
varnames_J <- c()
for(i in seq(nrow(upt_J))){
	for(j in seq(ncol(upt_J))){
		if(upt_J[i,j]){
			varnames_J <- c(varnames_J, paste0("Var", i, "-Var", j))
		}
	}
}

# Matrix => Vector
J_cont <- J_cont[upt_J]
J_dapt <- J_dapt[upt_J]

# Re-ordering
J_diff <- J_dapt - J_cont
order_J_diff <- order(unlist(J_diff), decreasing=TRUE)
J_cont <- J_cont[order_J_diff]
J_dapt <- J_dapt[order_J_diff]
J_diff <- J_diff[order_J_diff]
varnames_J <- varnames_J[order_J_diff]

# Data for ggplot2
gdata <- data.frame(value=c(J_cont, J_dapt),
        group=c(rep("cont", length=length(J_cont)),
            rep("dapt", length=length(J_dapt))),
        variable=factor(varnames_J, level=varnames_J))

# Plot
g <- ggplot() + geom_line(aes(x=variable, y=value, colour=group, group=group), linewidth=1, data = gdata, stat="identity") + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(outfile, plot=g, width=12, height=4.5)
