##### Paper plots: Example 9

library(igraph)

## 1. fully connected graph
G = cbind(c(0,2,2,2), c(2,0,2,2), c(2,2,0,2), c(2,2,2,0))
Gamma2Graph(G)

## 2. star-shaped graph
G = cbind(c(0,1,1,1), c(1,0,2,2), c(1,2,0,2), c(1,2,2,0))
Gamma2Graph(G)

## 3. Non-decomposable graph
G = cbind(c(0,1.5,1.5,2), c(1.5,0,2,1.5), c(1.5,2,0,1.5), c(2,1.5,1.5,0))
Gamma2Sigma(G, full=TRUE)
Gamma2Graph(G)



d=4
par=G
no.simu=8
model="HR"

if (model=="HR") {
  stopifnot(is.matrix(par))
  Gamma = par
  stopifnot(nrow(Gamma) == d & ncol(Gamma) == d)
  cov.mat <- Gamma2Sigma(Gamma, k=1, full=FALSE)
  chol.mat <- matrix(0,d,d)
  chol.mat[-1,-1] <- chol(cov.mat) ## add warning if cannot chol()
}


### Internal: simulates HR extremal functions on a tree
simu_px_tree_HR <- function(no.simu=1, G.vec, A) {
  res <- exp(A %*% matrix(rnorm(length(G.vec)*no.simu, mean= -G.vec/2, sd=sqrt(G.vec)), ncol=no.simu))
  return(t(res))
}

### Internal: simulates logistic extremal functions on a tree
simu_px_tree_logistic <- function(no.simu=1, idx, nb.edges, theta, A) {
  stopifnot(length(idx)==1 || length(idx)==no.simu)
  res       <- exp(A[[idx]] %*% log(matrix(1/gamma(1-theta)*(-log(runif(no.simu*nb.edges)))^(-theta) /
                                             (1/gamma(1-theta)*rgamma(no.simu*nb.edges,shape=1-theta)^(-theta)), ncol=no.simu)))
  return(t(res))
}

### Internal: simulates Dirichlet extremal functions on a tree
simu_px_tree_dirichlet <- function(no.simu=1, alpha.start, alpha.end, A) {
  e = length(alpha.start)
  shape.start = matrix(alpha.start + 1, nrow=e, ncol=no.simu)
  rate.start = matrix(alpha.start, nrow=e, ncol=no.simu)
  shape.end = matrix(alpha.end, nrow=e, ncol=no.simu)
  rate.end = matrix(alpha.end, nrow=e, ncol=no.simu)
  sim.start = matrix(rgamma(e*no.simu, shape=shape.start, rate=rate.start), nrow=e, ncol=no.simu)
  sim.end = matrix(rgamma(e*no.simu, shape=shape.end, rate=rate.end), nrow=e, ncol=no.simu)
  res       <- exp(A %*% log(sim.end / sim.start))
  return(t(res))
}


