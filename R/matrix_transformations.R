
#' Transformation of matrix to graph object
#'
#' Transforms a \eGamma or \eTheta matrix to an [`igraph::graph`] object for
#' the corresponding Huesler--Reiss extremal graphical model.
#'
#' @param Gamma Numeric \dxd variogram matrix.
#' @param tol Numeric scalar, entries in the precision matrix with absolute value
#' smaller than this are considered to be zero.
#'
#' @return Graph object from `igraph` package. An undirected graph.
#'
#' @examples
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' Gamma2graph(Gamma)
#'
#' @family MatrixTransformations
#' @seealso [get_large_tol()]
#' @rdname Gamma2graph
#' @export
Gamma2graph <- function(Gamma, tol=get_large_tol()){
  Theta2graph(Gamma2Theta(Gamma), tol=tol)
}
#' @param Theta Numeric \dxd precision matrix.
#' @rdname Gamma2graph
#' @export
Theta2graph <- function(Theta, tol=get_large_tol()){
  A <- 1*(abs(Theta) > tol)
  graph <- igraph::graph_from_adjacency_matrix(
    A,
    mode = "undirected",
    diag = FALSE
  )
  return(graph)
}

#' Transformation of a partial matrix to a graph
#' 
#' Creates a graph that has edges in entries corresponding to non-NA entries
#' in Gamma.
#' 
#' @param Matrix A matrix with NA entries
#' 
#' @return An `igraph::graph` object
#' @keywords internal
partialMatrixToGraph <- function(Matrix){
  A <- !is.na(Matrix)
  graph <- igraph::graph_from_adjacency_matrix(
    A,
    mode='undirected',
    diag = FALSE
  )
  return(graph)
}



#' Transformation of \eSigma and \eSigmaK matrix to \eGamma matrix
#'
#' Transforms the \eSigmaK matrix from the definition of a
#' Huesler--Reiss distribution to the corresponding \eGamma matrix.
#'
#' @param Sigma Numeric \dxd covariance matrix \eSigma, or \d1xd1 covariance
#' matrix \eSigmaK from the definition of a Huesler--Reiss distribution.
#' @param k `NULL` or integer between `1` and `d`.
#' If `NULL`, `Sigma` is interpreted as \eSigma, else
#' indicates which matrix \eSigmaK is given.
#'
#' @return Numeric \dxd \eGamma matrix.
#'
#' @examples
#' Sigma1 <- rbind(
#'   c(1.5, 0.5, 1),
#'   c(0.5, 1.5, 1),
#'   c(1, 1, 2)
#' )
#' Sigma2Gamma(Sigma1, k = 1, full = FALSE)
#' 
#' @family MatrixTransformations
#' @export
Sigma2Gamma <- function(Sigma, k = NULL) {
  if(is.null(k)) {
    # Sigma is already d x d
    S_full <- Sigma
  } else {
    # S is (d-1) x (d-1) => fill up with zeros
    d <- NROW(Sigma) + 1
    S_full <- matrix(0, d, d)
    S_full[-k ,-k] <- Sigma
  }

  # compute Gamma
  One <- rep(1, times = ncol(S_full))
  D <- diag(S_full)
  Gamma <- One %*% t(D) + D %*% t(One) - 2 * S_full

  return(Gamma)
}



#' Transformation of \eGamma matrix to \eSigma or \eSigmaK matrix
#'
#' Transforms the `Gamma` matrix from the definition of a
#' Huesler--Reiss distribution to the corresponding
#' \eSigma or \eSigmaK matrix.
#'
#' @param Gamma Numeric \dxd variogram matrix.
#' @param k `NULL` (default) or an integer between `1` and `d`.
#' Indicates which matrix \eSigma, or \eSigmaK
#' should be produced.
#' @param full Logical. If true, then the `k`th row and column in \eSigmaK
#' are included and the function returns a \dxd matrix.
#' By default, `full = FALSE`.
#'
#' @details
#' Every \dxd `Gamma` matrix in the definition of a
#' Huesler--Reiss distribution can be transformed into a
#' \d1xd1 \eSigmaK matrix,
#' for any `k` from `1` to `d`. The inverse of \eSigmaK
#' contains the graph structure corresponding to the Huesler--Reiss distribution
#' with parameter matrix `Gamma`. If `full = TRUE`, then \eSigmaK
#' is returned as a \dxd matrix with additional `k`th row and column
#' that contain zeros.
#' For details see \insertCite{eng2019;textual}{graphicalExtremes} and
#' \insertCite{hen2022;textual}{graphicalExtremes}.
#'
#' @return Numeric \eSigmaK matrix of size \d1xd1 if
#' `full = FALSE`, and \eSigma of size \dxd if `full = TRUE`.
#'
#' @examples
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#' Gamma2Sigma(Gamma, k = 1, full = FALSE)
#' 
#' @references \insertAllCited{}
#' 
#' @family MatrixTransformations
#' @export
Gamma2Sigma <- function(Gamma, k = NULL, full = FALSE) {
  d <- ncol(Gamma)
  if(is.null(k)) {
    ID <- diag(d) - matrix(1/d, d, d)
    ID %*% (-1/2 * Gamma) %*% ID
  } else if (full) {
    1 / 2 * (matrix(rep(Gamma[, k], d), ncol = d, nrow = d) +
      t(matrix(rep(Gamma[, k], d), ncol = d, nrow = d)) - Gamma)
  } else {
    1 / 2 * (matrix(rep(Gamma[-k, k], d - 1), ncol = d - 1, nrow = d - 1) +
      t(matrix(rep(Gamma[-k, k], d - 1), ncol = d - 1, nrow = d - 1)) - Gamma[-k, -k])
  }
}


#' Transformation of \eGamma matrix to \eTheta matrix
#'
#' Transforms the variogram matrix (\eGamma) from the definition of a
#' Huesler--Reiss distribution to the corresponding precision matrix
#' (\eTheta or \eThetaK).
#'
#'
#' @param Gamma Numeric \dxd variogram matrix.
#' @param k `NULL` or integer between 1 and d. If this is `NULL`, the
#' \dxd matrix \eTheta is produced, otherwise
#' the specified \eqn{(d-1) \times (d-1)}{(d-1) x (d-1)} matrix \eThetaK.
#' 
#' @details
#' Every \dxd `Gamma` matrix in the definition of a
#' Huesler--Reiss distribution can be transformed into a
#' \dxd \eTheta matrix, which
#' contains the graph structure corresponding to the Huesler--Reiss distribution
#' with parameter matrix `Gamma`.
#'
#' @return Numeric \eSigmaK matrix of size \d1xd1 if
#' `full = FALSE`, and of size \dxd if `full = TRUE`.
#'
#' @examples
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#' Gamma2Theta(Gamma)
#' 
#' @references \insertAllCited{}
#' 
#' @family MatrixTransformations
#' @export
Gamma2Theta <- function(Gamma, k=NULL) {
  d <- ncol(Gamma)
  if(is.null(k)){
    ID <- diag(d) - matrix(1/d, d, d)
    S <- ID %*% (-1/2 * Gamma) %*% ID
    Theta <- corpcor::pseudoinverse(S)
    return(Theta)
  }
  Sigma_k <- Gamma2Sigma(Gamma, k)
  Theta_k <- chol2inv(chol(Sigma_k))
  return(Theta_k)
}

#' Transformation of \eGamma matrix to \eTheta matrix
#' 
#' Transforms a precision matrix (\eTheta or \eThetaK)
#' to the corresponding variogram matrix.
#' 
#' @param Theta Numeric \dxd matrix (if `k` is `NULL`)
#' or \eqn{(d-1) \times (d-1)}{(d-1) x (d-1)} matrix (if `k` is a number).
#' @param k `NULL` or integer between 1 and d.
#' If this is `NULL` the input `Theta` is interpreted as a \dxd
#' precision matrix \eTheta, otherwise as \eThetaK.
#' 
#' @return The \dxd variogram matrix implied by `Theta`.
#' 
#' @family MatrixTransformations
#' @export
Theta2Gamma <- function(Theta, k=NULL) {
  if(is.null(k)){
    Sigma <- corpcor::pseudoinverse(Theta)
    return(Sigma2Gamma(Sigma))
  }
  Sigma_k <- chol2inv(chol(Theta))
  return(Sigma2Gamma(Sigma_k, k))
}


#' Create Gamma or Theta from vector
#'
#' This function takes the parameters in the vector `par`
#' (upper triangular Gamma/Theta matrix) and returns the full Gamma/Theta.
#'
#' @param par Numeric vector with `d` elements.
#' Upper triangular part of a Gamma/Theta matrix.
#' @param allowMatrix If `TRUE` and `par` is already a matrix, return it as is.
#' @param allowNull If `TRUE` and `par` is NULL, return NULL.
#' @param zeroRowSums If `TRUE` the diagonal is set to (-1) times the rowSums.
#'
#' @return Numeric matrix \dxd. Full Gamma/Theta matrix corresponding to `par`.
#' 
#' @rdname par2Matrix
#' @keywords internal
par2Matrix <- function(par, allowMatrix = FALSE, allowNull = FALSE, zeroRowSums = FALSE){
  # Check for forbidden/trivial input
  if(is.matrix(par)){
    if(allowMatrix){
      return(par)
    }
    stop('`par` must not be a matrix (unless allowMatrix=TRUE).')
  }
  if(is.null(par)){
    if(allowNull){
      return(NULL)
    }
    stop('`par` must not be NULL (unless allowNull=TRUE).')
  }
  # Compute dimension of matrix
  d <- round(1 / 2 + sqrt(1 / 4 + 2 * length(par)))
  if (d*(d-1)/2 != length(par)) {
    stop("The length of par does not agree with the upper triangle of any square matrix.")
  }
  # Create matrix
  M <- matrix(0, nrow = d, ncol = d)
  M[upper.tri(M)] <- par
  M <- M + t(M)
  if(zeroRowSums){
    diag(M) <- (-1) * rowSums(M)
  }
  return(M)  
}
#' @rdname par2Matrix
par2Gamma <- function(par, allowMatrix = FALSE, allowNull = FALSE){
  par2Matrix(par, allowMatrix, allowNull, zeroRowSums = FALSE)
}
#' @rdname par2Matrix
par2Theta <- function(par, allowMatrix = FALSE, allowNull = FALSE){
  par2Matrix(par, allowMatrix, allowNull, zeroRowSums = TRUE)
}
#' @rdname par2Matrix
#' @param allowVector If `TRUE` and `M` is already a vector, return it as is.
#' @param M Matrix
#' 
#' @return Upper triangular part of `M` (or `M` itself/NULL if allowed)
matrix2par <- function(M, allowVector = FALSE, allowNull = FALSE){
  if(is.null(M)){
    if(allowNULL){
      return(NULL)
    }
    stop('Matrix `M` must not be NULL (unless allowNull=TRUE).')
  }
  if(allowVector && is.vector(M)){
    return(M)
  }
  return(upper.tri.val(M))
}


#' Extract upper triangular part of \eGamma
#'
#' This function returns a vector containing the upper triangular part
#' of the matrix `Gamma`. If `Gamma` is already a vector, it returns
#' it as it is.
#'
#' @param Gamma Numeric \dxd variogram matrix.
#'
#' @return Numeric vector with `d` elements.
#' The upper triangular part of the given `Gamma` matrix.
#'
#' @keywords internal
Gamma2par <- function(Gamma) {
  if(is.matrix(Gamma)) {
    return(upper.tri.val(Gamma))
  }
  return(Gamma)
}



#' Transformation of extremal correlation \eChi to the Huesler--Reiss variogram \eGamma
#'
#' Transforms the extremal correlation \eChi into the `Gamma` matrix
#' from the definition of a Huesler--Reiss distribution.
#'
#' @param chi Numeric or matrix, with entries
#' between 0 and 1.
#'
#' @return Numeric or matrix. The \eGamma parameters in the Huesler--Reiss
#' distribution.
#'
#' @details
#' The formula for transformation from `chi` to \eGamma that is applied element-wise is
#' \deqn{\Gamma = (2 \Phi^{-1}(1 - 0.5 \chi))^2,}
#' where \eqn{\Phi^{-1}} is the inverse of the standard normal distribution function.
#' This is the inverse of [Gamma2chi()].
#'
#' @export
chi2Gamma <- function(chi) {
  if (any(is_less(chi, 0) | is_greater(chi, 1))) {
    stop("The argument chi must be between 0 and 1.")
  }
  Gamma <- (2 * stats::qnorm(1 - 0.5 * chi))^2
  return(Gamma)
}


#' Transformation of the Huesler--Reiss variogram \eGamma to extremal correlation \eChi
#'
#' Transforms the `Gamma` matrix from the definition of a Huesler--Reiss
#' distribution into the corresponding extremal correlation \eChi.
#'
#' @details
#' The formula for transformation from `Gamma` to \eChi that is applied element-wise is
#' \deqn{\chi = 2 - 2 \Phi(\sqrt{\Gamma} / 2),}{\chi = 2 - 2 \Phi(sqrt(\Gamma) / 2),}
#' where \eqn{\Phi} is the standard normal distribution function.
#' This is the inverse of [chi2Gamma()].
#'
#'
#' @param Gamma Numeric or matrix, with positive entries.
#'
#' @return Numeric or matrix. The extremal correlation coefficient.
#'
#' @export
Gamma2chi <- function(Gamma) {
  chi <- 2 - 2 * stats::pnorm(sqrt(Gamma) / 2)
  return(chi)
}


#' Compute theoretical \eChi in 3D
#'
#' Computes the theoretical \eChi coefficient in 3 dimensions.
#'
#' @param Gamma Numeric \eXTimesY{3}{3} matrix.
#'
#' @return The 3-dimensional \eChi coefficient, i.e.,
#' the extremal correlation coefficient for the HR distribution. Note that
#' \eqn{0 \leq \chi \leq 1}.
#'
#' @keywords internal
Gamma2chi_3D <- function(Gamma) {
  d <- dim_Gamma(Gamma)

  if (d != 3) {
    stop("Gamma must be a 3 x 3 matrix.")
  }
  res <- 3 - V_HR(x = rep(1, times = 2), Gamma = Gamma[c(1, 2), c(1, 2)]) -
    V_HR(x = rep(1, times = 2), Gamma = Gamma[c(1, 3), c(1, 3)]) -
    V_HR(x = rep(1, times = 2), Gamma = Gamma[c(2, 3), c(2, 3)]) +
    V_HR(x = rep(1, times = 3), Gamma = Gamma)
  return(res)
}

