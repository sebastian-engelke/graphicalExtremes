

nA <- 2
nC <- 3
nB <- 4
nV <- nA + nB + nC

A <- 1:nA
AC <- 1:(nA + nC)
C <- (nA + 1):(nA + nC)
BC <- (nA + 1):nV
B <- (nA + nC + 1):nV
V <- 1:nV

G0 <- generate_random_Gamma(nV)

G <- matrix(NA, nV, nV)
G[AC,AC] <- G0[AC,AC]
G[BC,BC] <- G0[BC,BC]

G1 <- graphicalExtremes::complete_Gamma(G)
Gnew1 <- G1[A, B]

make_eD <- function(V, D){
  eD <- numeric(length(V))
  eD[D] <- 1 / length(D)
  return(eD)
}
make_Pv <- function(v){
  oneVec <- numeric(length(v)) + 1
  return(diag(length(v)) - oneVec %*% t(v))
}
make_PD <- function(V, D){
  make_Pv(make_eD(V, D))
}


makeSD <- function(G, D){
  G2 <- G
  G[is.na(G)] <- 1
  G2[is.na(G2)] <- 10
  eD <- make_eD(V, D)
  PD <- make_Pv(eD)
  S <- PD %*% (-G/2) %*% t(PD)
  S2 <- PD %*% (-G2/2) %*% t(PD)
  diffInd <- (abs(S - S2)) > 1e-6
  S[diffInd] <- NA
  return(S)
}

oneVec <- numeric(nV) + 1
eD <- make_eD(V, C)
PD <- make_Pv(eD)

# S <- PD %*% (-G0/2) %*% t(PD)

SD <- makeSD(G, C)

S1 <- makeSD(G1, C)

ThetaC <- corpcor::pseudoinverse(SD[C, C, drop=FALSE])
SDnew <- SD[A, C] %*% ThetaC %*% SD[C, B]

S2 <- SD
S2[A,B] <- SDnew
S2[B,A] <- t(SDnew)
G2 <- Sigma2Gamma(S2)


