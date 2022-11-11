

#' Fit value(s) in interval
#' 
#' Fit value(s) in interval, all arguments are recycled where necessary.
#' 
#' @param x Numeric vector
#' @param xMin Numeric vector
#' @param xMax Numeric vector
#' @return Numeric vector
fitInInterval <- function(x, xMin=-Inf, xMax=Inf){
  if(any(xMax<xMin)){
    stop('Make sure that xMax>=xMin!')
  }
  x <- pmax(x, xMin)
  x <- pmin(x, xMax)
  return(x)
}


# Replaces a principal submatrix of a Gamma matrix, preserving definiteness.
# Other entries are kept "heuristically similar".
# Zeros (instead of NAs) on the diagonal of G.fix indicate which principal submatrix to replace.
# If not the entire submatrix in G.fix is specified, [complete_Gamma_decomposable] is used
replaceGammaSubMatrix <- function(G.est, G.fix){
  # check which are fixed
  ind <- which(!is.na(diag(G.fix)))
  if(length(ind)==0){
    return(G.est)
  }
  if(any(is.na(G.fix[ind, ind]))){
    G.fix[ind, ind] <- complete_Gamma(G.fix[ind, ind], allowed_graph_type = 'decomposable')
  }
  if(length(ind) == ncol(G.est)){
    return(G.fix)
  }
  
  # naive attempt:
  G1 <- G.est
  G1[ind, ind] <- G.fix[ind, ind]
  if(is_sym_cnd(G1)){
    return(G1)
  }

  # heuristic, but safe solution:
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
    warning('Matrix not symmetric (up to tolerance: ', tol, ')!')
  }
  (M + t(M))/2
}

is_symmetric <- function(M, tol=1e-12){
  max(M - t(M)) < tol
}

is_cnd <- function(M, tol=1e-12){
  d <- nrow(M)
  if(d == 0){
    return(TRUE)
  }
  if(d == 1){
    return(abs(M[1,1]) < tol)
  }
  if(any(upper.tri.val(M) <= 0)){
    return(FALSE)
  }
  Sk <- Gamma2Sigma(M, k=1)
  eig <- eigen(Sk, symmetric = TRUE, only.values = TRUE)$values
  return(abs(eig[d-1]) > 0)
}

is_sym_cnd <- function(M, tol=1e-12){
  is_symmetric(M, tol) && is_cnd(M, tol)
}

is_valid_Theta <- function(M, tol=1e-9){
  if(!is_symmetric(M, tol)){
    return(FALSE)
  }
  d <- nrow(M)
  if(d == 0){
    return(TRUE)
  }
  if(d == 1){
    return(abs(M[1,1]) < tol)
  }
  if(any(abs(rowSums(M)) > tol)){
    return(FALSE)
  }
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if(abs(eig[d]) > tol || abs(eig[d-1]) < tol){
    return(FALSE)
  }
  return(TRUE)
}


upper.tri.val <- function(M, diag=FALSE){
  M[upper.tri(M, diag)]
}


DEFAULT_TOL <- .Machine$double.eps^0.5
is_eq <- function(a, b, tol=NULL) {
  if(is.null(tol)){
    tol <- DEFAULT_TOL
  }
  abs(a - b) < tol
}
is_greater <- function(a, b, tol=NULL) {
  if(is.null(tol)){
    tol <- DEFAULT_TOL
  }
  a - b > tol
}
is_less <- function(a, b, tol=NULL) {
  is_greater(b, a, tol)
}
is_leq <- function(a, b, tol=NULL) {
  !is_greater(a, b, tol)
}
is_geq <- function(a, b, tol=NULL) {
  !is_less(a, b, tol)
}



