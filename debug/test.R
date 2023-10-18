
# devtools::load_all()

options(warn=1)

eps <- 1e-15
eps2 <- 1e-13
eps3 <- 0.0004

Mij <- function(i, j, d){
  M <- matrix(0, d, d)
  M[i,j] <- 1
  return(M)
}

emptyMatrix <- matrix(0, 0, 0)

oneMatrix0 <- matrix(0, 1, 1)
oneMatrix1 <- matrix(1, 1, 1)
oneMatrixEps <- matrix(eps, 1, 1)

G1 <- matrix(c(0, 1, 1, 0), 2, 2)
G2 <- G1 + diag(2) * eps
G3 <- G1 + Mij(1, 2, 2) * eps2 + diag(2) * eps2

G3c <- checkGamma(G3, tol = 0, alert=TRUE)

G3c <- checkGamma(G3c, tol = 0, alert=TRUE)
