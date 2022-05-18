
set.seed(1)

d <- 20

G <- generate_random_model(d)$Gamma

g <- generate_random_connected_graph(d, m = d*4)

ret <- complete_Gamma_general(G, g, N=10000, saveDetails=TRUE)
GammaList <- ret$GammaList
G2 <- GammaList[[length(GammaList)]]

P2 <- Gamma2Theta(G2)
A <- igraph::as_adjacency_matrix(g, sparse = FALSE)
B <- igraph::as_adjacency_matrix(igraph::complementer(g), sparse = FALSE)
B <- B * upper.tri(B)

m <- sum(B)

print(max(abs(P2[B == 1])))
print(max(abs(G2 - G)[A == 1]))

ThetaList <- lapply(GammaList, Gamma2Theta)

eVec <- c()
iVec <- c()
pMat <- matrix(0, m, 0)
for(i in seq_along(ThetaList)){
  P3 <- ThetaList[[i]]
  eVec[i] <- max(abs(P3 * B))
  iVec[i] <- which.max(P3 * B)
  pMat <- cbind(pMat, P3[B == 1])
}

x <- seq_along(eVec)

nOrd <- 4

ord <- order(tabulate(iVec), decreasing = TRUE)

# plot(cbind(x, x, x), log(cbind(eVec, p1Vec, p2Vec)), type='l')

mx <- log(max(abs(c(eVec, pMat))))
mn <- log(min(abs(c(eVec, pMat))))

plot(x, log(eVec), type='l', ylim=c(mn, mx), lwd = 2)

for(i in seq_len(m)){
  lines(x, log(abs(pMat[i,])), type='l', col=i+1)
}

lines(x, log(eVec), type='l', col = 1, lwd = 2)





