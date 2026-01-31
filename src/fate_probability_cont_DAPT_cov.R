source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
infile3 <- commandArgs(trailingOnly=TRUE)[3]
infile4 <- commandArgs(trailingOnly=TRUE)[4]
infile5 <- commandArgs(trailingOnly=TRUE)[5]
infile6 <- commandArgs(trailingOnly=TRUE)[6]
outfile1 <- commandArgs(trailingOnly=TRUE)[7]
outfile2 <- commandArgs(trailingOnly=TRUE)[8]
outfile3 <- commandArgs(trailingOnly=TRUE)[9]
outfile4 <- commandArgs(trailingOnly=TRUE)[10]
# infile1 <- 'plot/hpbase/cont_cov/Landscaper/Basin.tsv'
# infile2 <- 'plot/hpbase/DAPT_cov/Landscaper/Basin.tsv'
# infile3 <- 'plot/hpbase/cont_cov/P_metropolis.tsv'
# infile4 <- 'plot/hpbase/DAPT_cov/P_metropolis.tsv'
# infile5 <- 'plot/hpbase/cont_cov/P_glauber.tsv'
# infile6 <- 'plot/hpbase/DAPT_cov/P_glauber.tsv'

# Loading
Basin_cont <- as.vector(unlist(read.table(infile1, header=FALSE)))
Basin_DAPT <- as.vector(unlist(read.table(infile2, header=FALSE)))
Basin <- union(Basin_cont, Basin_DAPT)

# cont
P_m <- as.matrix(read.table(infile3, header=FALSE))
P_g <- as.matrix(read.table(infile5, header=FALSE))

## Metropolis
res_m <- absorption_probabilities(P_m, absorbing = Basin, reg = 1e-4)
H_m <- fate_entropy(res_m$F, base = 2)
argmax_m <- fate_argmax(res_m$F)

## Glauber
res_g <- absorption_probabilities(P_g, absorbing = Basin, reg = 1e-4)
H_g <- fate_entropy(res_g$F, base = 2)
argmax_g <- fate_argmax(res_g$F)

## Save
save(res_m, H_m, argmax_m, file=outfile1)
save(res_g, H_g, argmax_g, file=outfile2)


# DAPT
P_m <- as.matrix(read.table(infile4, header=FALSE))
P_g <- as.matrix(read.table(infile5, header=FALSE))

## Metropolis
res_m <- absorption_probabilities(P_m, absorbing = Basin, reg = 1e-4)
H_m <- fate_entropy(res_m$F, base = 2)
argmax_m <- fate_argmax(res_m$F)

## Glauber
res_g <- absorption_probabilities(P_g, absorbing = Basin, reg = 1e-4)
H_g <- fate_entropy(res_g$F, base = 2)
argmax_g <- fate_argmax(res_g$F)

## Save
save(res_m, H_m, argmax_m, file=outfile3)
save(res_g, H_g, argmax_g, file=outfile4)
