source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
outfile1 <- commandArgs(trailingOnly=TRUE)[4]
outfile2 <- commandArgs(trailingOnly=TRUE)[5]
# infile1 <- 'plot/hpbase/integrated/Landscaper/Basin.tsv'
# infile2 <- 'plot/hpbase/integrated/P_metropolis.tsv'
# infile3 <- 'plot/hpbase/integrated/P_glauber.tsv'

# Loading
Basin <- as.vector(unlist(read.table(infile1, header=FALSE)))
P_m <- as.matrix(read.table(infile2, header=FALSE))
P_g <- as.matrix(read.table(infile3, header=FALSE))

# Metropolis
res_m <- absorption_probabilities(P_m, absorbing = Basin, reg = 1e-4)
H_m <- fate_entropy(res_m$F, base = 2)
argmax_m <- fate_argmax(res_m$F)

# Glauber
res_g <- absorption_probabilities(P_g, absorbing = Basin, reg = 1e-4)
H_g <- fate_entropy(res_g$F, base = 2)
argmax_g <- fate_argmax(res_g$F)

# Save
save(res_m, H_m, argmax_m, file=outfile1)
save(res_g, H_g, argmax_g, file=outfile2)
