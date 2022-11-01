
set.seed(1)

d <- 10
N <- 100
nLines <- 10

g_comp <- igraph::make_full_graph(d)
G <- generate_random_graphical_Gamma(g_comp)

g <- generate_random_connected_graph(d, m = d*4)

ret <- complete_Gamma_general(G, g, N=N-1, saveDetails=TRUE, check_tol = N)

# This step requires a modified version of `graphicalExtremes`:
GammaList <- ret$GammaList

A <- (igraph::as_adjacency_matrix(g, sparse = FALSE) == 1)
B <- (!A) & upper.tri(A)

m <- sum(B)


ThetaList <- lapply(GammaList, Gamma2Theta)

eVec <- numeric(N)
iVec <- numeric(N)
pMat <- matrix(0, m, N)
for(i in seq_along(ThetaList)){
  PB <- ThetaList[[i]][B]
  pMat[,i] <- PB
  iVec[i] <- which.max(abs(PB))
  eVec[i] <- abs(PB[iVec[i]])
}

ord <- order(tabulate(iVec), decreasing = TRUE)

# plot(cbind(x, x, x), log(cbind(eVec, p1Vec, p2Vec)), type='l')

mx <- log10(max(abs(c(eVec, pMat))))
mn <- log10(min(abs(c(eVec, pMat))))

plot(NULL, xlab='Iteration', ylab='$\\log_{10}(\\Theta_{ij})$', xlim=c(0, N), ylim=c(mn, mx))

for(i in seq_len(min(m, nLines))){
  lines(log10(abs(pMat[ord[i],])), type='l', col=i+1)
}

lines(log10(eVec), type='l', col = 1, lwd = 2)





