##### Paper plots: Example 9

## 1. fully connected graph
G = cbind(c(0,2,2,2), c(2,0,2,2), c(2,2,0,2), c(2,2,2,0))
Gamma2Graph(G)

## 2. star-shaped graph
G = cbind(c(0,1,1,1), c(1,0,2,2), c(1,2,0,2), c(1,2,2,0))
Gamma2Graph(G)

## 3. Non-decomposable graph
G = cbind(c(0,1.5,1.5,2), c(1.5,0,2,1.5), c(1.5,2,0,1.5), c(2,1.5,1.5,0))
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
