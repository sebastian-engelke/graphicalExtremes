
d <- 5
oneVec <- matrix(1, d, 1)
eye <- diag(d)
P <- diag(d) - oneVec %*% t(oneVec) * (1/d)

m <- generate_random_model(d)

G <- m$Gamma

S <- Gamma2Sigma(G)

T <- corpcor::pseudoinverse(S)

D <- -1/2*G

ss <- solve(D, oneVec)
s <- ss / sum(ss)

F <- (eye - oneVec %*% t(s)) %*% D %*% (eye - s %*% t(oneVec))

makeX <- function(A){
  V <- eigen(A)$vector
  w <- eigen(A)$values
  w[abs(w) < 1e-8] <- 0
  V <- V %*% diag(sqrt(w))
  V[,abs(w) < 1e-8] <- 0
  return(V)
}

X <- makeX(F)

l <- diag(1/sqrt(diag(F)))

Y <- X %*% l

F2 <- Y %*% t(Y)

S2 <- P %*% F2 %*% P

G2 <- Sigma2Gamma(S2)

T2 <- Gamma2Theta(G2)


