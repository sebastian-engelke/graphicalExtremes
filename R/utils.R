fitInInterval <- function(x, xMin=-Inf, xMax=Inf){
  if(any(xMax<xMin)){
    stop('Make sure that xMax>=xMin!')
  }
  x <- pmax(x, xMin)
  x <- pmin(x, xMax)
  return(x)
}


# replaces part of a Gamma matrix, keeping definiteness
replaceGammaSubMatrix <- function(G.est, G.fix){
  ind <- which(!is.na(diag(G.fix)))
  if(length(ind)==0){
    return(G.est)
  }
  k <- ind[1]
  M.est <- Gamma2Sigma(G.est, k)
  M.fix <- Gamma2Sigma(G.fix, k)
  M <- replaceSpdSubmatrix(M.est, M.fix)
  G <- Sigma2Gamma(M, k=k)
  return(G)
}

# replaces part of a positive definite matrix, keeping definiteness
replaceSpdSubmatrix <- function(M.est, M.fix){
  indFix <- !is.na(diag(M.fix))
  indMod <- !indFix

  indFix <- which(indFix)
  indMod <- which(indMod)

  if(length(indFix) == 0){
    return(M.est)
  } else if(length(indMod) == 0){
    return(M.fix)
  }

  ind <- c(indFix, indMod)

  M2 <- M.est[ind, ind]

  indFix2 <- seq_along(indFix)
  indMod2 <- seq_along(indMod) + length(indFix2)

  C <- M.fix[indFix, indFix]

  L <- chol(M2)
  LC <- chol(C)

  L[indFix2, indFix2] <- LC

  M2 <- t(L) %*% L

  M <- matrix(0, NROW(M2), NROW(M2))

  M[ind, ind] <- M2

  return(M)
}

rdunif <- function(n, a, b){
  a + floor((b - a + 1) * runif(n))
}


ensure_symmetry <- function(M, tol=1e-6){
  if(max(abs(M - t(M))) > tol){
    warning(paste0(
      'Matrix not symmetric (up to tolerance: ',
      max(abs(M - t(M))),
      ' > ',
      tol,
      ')!'
    ))
  }
  (M + t(M))/2
}

is_sym_cnd <- function(M){
  if(!matrixcalc::is.symmetric.matrix(M)){
    return(FALSE)
  }
  Sk <- Gamma2Sigma(M, k=1)
  return(matrixcalc::is.positive.definite(Sk))
}
