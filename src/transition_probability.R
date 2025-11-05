source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
# infile1 <- 'plot/hpbase/integrated/Landscaper/E.tsv'
# infile2 <- 'plot/hpbase/integrated/Landscaper/Allstates.tsv'

# Loading
E <- unlist(read.table(infile1, header=FALSE))
Allstates <- as.matrix(read.table(infile2, header=FALSE))

# Transition Probability
P_metropolis <- transition_matrix_from_E(E, Allstates, beta = 1, kernel = "metropolis")
P_glauber <- transition_matrix_from_E(E, Allstates, beta = 1, kernel = "glauber")

# Output
write.table(P_metropolis, file=outfile1,
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
write.table(P_glauber, file=outfile2,
    sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)
