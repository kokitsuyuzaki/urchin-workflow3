source("src/Functions.R")

# Parameter
infile1 <- commandArgs(trailingOnly=TRUE)[1]
infile2 <- commandArgs(trailingOnly=TRUE)[2]
outfile1 <- commandArgs(trailingOnly=TRUE)[3]
outfile2 <- commandArgs(trailingOnly=TRUE)[4]
# infile1 <- 'plot/hpbase/integrated/P_metropolis.tsv'
# infile2 <- 'plot/hpbase/integrated/P_glauber.tsv'

# Loading
P_m <- t(as.matrix(read.table(infile1, header=FALSE)))
P_g <- t(as.matrix(read.table(infile2, header=FALSE)))

# Random walk (metropolis)
P_m <- P_m %*% diag(rep(1, length=nrow(P_m)))

iter <- 1
tol <- 100
while(iter <= 100 & tol > 1e-6){
    P_m_new <- P_m %*% P_m
    tol <- max(abs(P_m_new - P_m))
    P_m <- P_m_new
    iter <- iter + 1
}
print(iter)

# Random walk (glauber)
P_g <- P_g %*% diag(rep(1, length=nrow(P_g)))

iter <- 1
tol <- 100
while(iter <= 100 & tol > 1e-6){
    P_g_new <- P_g %*% P_g
    tol <- max(abs(P_g_new - P_g))
    P_g <- P_g_new
    iter <- iter + 1
}
print(iter)

# Save
write.table(P_m, file=outfile1,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
write.table(P_g, file=outfile2,
    quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
    