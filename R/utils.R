

#' Fit value(s) in interval
#' 
#' Fit value(s) in interval, all arguments are recycled where necessary.
#' 
#' @param x Numeric vector
#' @param xMin Numeric vector
#' @param xMax Numeric vector
#' @return Numeric vector
#' @keywords internal
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
# Zeros (instead of NAs) on the diagonal of `Gamma.fix` indicate which principal submatrix to replace.
# If not the entire submatrix in `Gamma.fix` is specified, [complete_Gamma_decomposable()] is used
replaceGammaSubMatrix <- function(Gamma.est, Gamma.fix){
  # check which are fixed
  ind <- which(!is.na(diag(Gamma.fix)))
  if(length(ind)==0){
    return(Gamma.est)
  }
  if(any(is.na(Gamma.fix[ind, ind]))){
    Gamma.fix[ind, ind] <- complete_Gamma(Gamma.fix[ind, ind], allowed_graph_type = 'decomposable')
  }
  if(length(ind) == ncol(Gamma.est)){
    return(Gamma.fix)
  }
  
  # naive attempt:
  G1 <- Gamma.est
  G1[ind, ind] <- Gamma.fix[ind, ind]
  if(is_sym_cnd(G1)){
    return(G1)
  }

  # heuristic, but safe solution:
  k <- ind[1]
  M.est <- Gamma2Sigma(Gamma.est, k, check = FALSE)
  M.fix <- Gamma2Sigma(Gamma.fix, k, check = FALSE)
  M <- replaceSpdSubmatrix(M.est, M.fix)
  G <- Sigma2Gamma(M, k=k, check = FALSE)
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
  a + floor((b - a + 1) * stats::runif(n))
}

#' Ensure numerical matrix symmetry/zero values
#' 
#' Ensures the symmetry of a square matrix by averaging it with its transpose.
#' Optionally verifies that the matrix was close to symmetric before.
#' 
#' @param M Numeric square matrix.
#' @param checkTol Positive scalar. If the maximum absolute difference between `M`
#' and `t(M)` is larger, show a warning.
#' 
#' @return The adjusted value of `M`.
#' 
#' @rdname ensure_matrix_symmetry_and_truncate_zeros
#' @export
ensure_matrix_symmetry <- function(M, checkTol=Inf){
  if(checkTol < Inf && max(abs(M - t(M))) > checkTol){
    warning('Matrix not symmetric (up to tolerance: ', checkTol, ')!')
  }
  (M + t(M))/2
}
#' @description Makes sure zeros are "numerically zero", by truncating all small values.
#' @param tol All entries with absolute value below this value are truncated to zero.
#' 
#' @rdname ensure_matrix_symmetry_and_truncate_zeros
#' @export
truncate_zeros <- function(M, tol=get_small_tol()){
  M[abs(M) < tol] <- 0
  M
}
#' @rdname ensure_matrix_symmetry_and_truncate_zeros
#' @export
ensure_matrix_symmetry_and_truncate_zeros <- function(M, tol=get_small_tol(), checkTol=Inf){
  M <- ensure_matrix_symmetry(M, checkTol)
  truncate_zeros(M, tol)
}


is_square <- function(M){
  is.matrix(M) && ncol(M) == nrow(M)
}

# Check if a matrix is symmetric up to tolerance, allowing NAs
is_symmetric <- function(M, tol=get_small_tol()){
  (
    is_square(M)
    && all(is.na(M) == t(is.na(M)))
    && max(M - t(M), na.rm = TRUE) < tol
  )
}

is_sym_cnd <- function(M, tol=get_small_tol()){
  # M must be symmetric
  if(!is_symmetric(M, tol)){
    return(FALSE)
  }
  d <- nrow(M)
  # Empty matrix is cnd
  if(d == 0){
    return(TRUE)
  }
  # Diagonal must be zeros
  if(any(abs(diag(M)) > tol)){
    return(FALSE)
  }
  # 1x1 matrix is just 0 => ok
  if(d == 1){
    return(TRUE)
  }
  # All entries must be positive
  if(any(upper.tri.val(M) <= 0)){
    return(FALSE)
  }
  # Check that Gamma2Sigma yields a pos. def. matrix
  # Is there a more elegant way to do this?
  Sk <- Gamma2Sigma(M, k=1, check=FALSE)
  eig <- eigen(Sk, symmetric = TRUE, only.values = TRUE)$values
  return(eig[d-1] > 0)
}

is_sym_pos_def <- function(M, tol=get_small_tol()){
  if(!is_symmetric(M, tol)){
    return(FALSE)
  }
  # If M is symmetric, pos.def. is equivalent to all positive eigenvalues
  d <- nrow(M)
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  return(eig[d] > 0)
}

is_valid_Theta <- function(M, tol=get_small_tol()){
  # Must be symmetric
  if(!is_symmetric(M, tol)){
    return(FALSE)
  }
  d <- nrow(M)
  # Empty matrix is valid
  if(d == 0){
    return(TRUE)
  }
  # Check that rowsums are zero
  if(any(abs(rowSums(M)) > tol)){
    return(FALSE)
  }
  # 1x1 matrix is just 0 => ok
  if(d == 1){
    return(TRUE)
  }
  # Check that there are d-1 positive and one zero eigenvalue
  eig <- eigen(M, symmetric = TRUE, only.values = TRUE)$values
  if(abs(eig[d]) > tol || abs(eig[d-1]) <= 0){
    return(FALSE)
  }
  return(TRUE)
}

pdet <- function(M, tol=get_small_tol()){
  ev <- eigen(M, only.values = TRUE)$values
  prod(ev[abs(ev) > tol])
}


upper.tri.val <- function(M, diag=FALSE){
  M[upper.tri(M, diag)]
}


is_eq <- function(a, b, tol=get_small_tol()) {
  abs(a - b) < tol
}
is_greater <- function(a, b, tol=get_small_tol()) {
  a - b > tol
}
is_less <- function(a, b, tol=get_small_tol()) {
  is_greater(b, a, tol)
}
is_leq <- function(a, b, tol=get_small_tol()) {
  !is_greater(a, b, tol)
}
is_geq <- function(a, b, tol=get_small_tol()) {
  !is_less(a, b, tol)
}


makeUnitVector <- function(d, k){
    ek <- numeric(d)
    ek[k] <- 1
    return(ek)
}
makeOneVec <- function(d){
    v <- rep(1, d)
}
makeProjK <- function(d, k){
    ek <- makeUnitVector(d, k)
    oneVec <- rep(1, d)
    diag(d) - ek %*% t(oneVec)
}



#' Convert indices to numerical indices
#' 
#' Converts (possibly) logical indices to numerical ones.
#' Also ensures unique indices and sorts them if specified.
#' 
#' @param ind The numerical or logical index vector
#' @param n Max numerical index (used if `ind` is logical and might be recycled)
#' @param unique Whether to keep every (numerical) index at most once
#' @param sort Whether to sort the numerical indices
#' 
#' @return A numerical index vector
#' 
#' @keywords internal
make_numeric_indices <- function(ind, n=NULL, unique=TRUE, sort=TRUE){
  if(is.logical(ind)){
    if(is.null(n)){
      n <- length(ind)
    }
    ind0 <- seq_len(n)
    ind <- ind0[ind]
  }
  ind <- as.integer(ind)
  if(unique){
    ind <- unique(ind)
  }
  if(sort){
    ind <- sort(ind)
  }
  return(ind)
}

