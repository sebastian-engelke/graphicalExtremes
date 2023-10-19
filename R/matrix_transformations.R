

#' @param Sigma Numeric \dxd or \d1xd1 covariance matrix.
#' @param Gamma Numeric \dxd variogram matrix.
#' @param Theta Numeric \dxd or \d1xd1 precision matrix.
#' @param k `NULL` if the input/output matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param k1 `NULL` if the input matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param k2 `NULL` if the output matrix is \eSigma/\eTheta.
#' Else, an integer between 1 and d indicating the value of k in \eSigmaK, \eThetaK.
#' @param full Logical. If `TRUE` and `!is.null(k)`,
#' the input/output matrix is a \dxd matrix with the kth row filled with zeros.
#' @param full1 Logical. If `TRUE` and `!is.null(k1)`,
#' the input is a \dxd matrix with the kth row filled with zeros.
#' @param full2 Logical. If `TRUE` and `!is.null(k2)`,
#' the output is a \dxd matrix with the kth row filled with zeros.
#' @param check Whether to check the inputs and call `ensure_matrix_symmetry_and_truncate_zeros`
#' on the outputs.
#' 
#' @name sharedParamsMatrixTransformations
NULL

#' Convert matrix to graph
#' 
#' Creates a graph object representing the graph structure implied by a parameter matrix.
#' 
#' @param M Partial matrix with `NA` entries indicating missing edges.
#' @param tol Numeric scalar. Entries in the precision matrix with absolute value
#' smaller than this are considered to be zero.
#' @inheritParams sharedParamsMatrixTransformations
#' 
#' @return An [`igraph::graph`] object.
#' 
#' @family matrixTransformations
#' @rdname matrix2graph
#' @export
Gamma2graph <- function(Gamma, tol=get_large_tol(), check=TRUE){
  Theta <- Gamma2Theta(Gamma, check=check)
  Theta2graph(Theta, tol=tol, check=FALSE)
}
#' @rdname matrix2graph
#' @export
Sigma2graph <- function(Sigma, tol=get_large_tol(), k=NULL, full=FALSE, check=TRUE){
  Theta <- Sigma2Theta(Sigma, k1=k, full1=full, check=check)
  Theta2graph(Theta)
}
#' @rdname matrix2graph
#' @export
Theta2graph <- function(Theta, tol=get_large_tol(), k=NULL, full=FALSE, check=TRUE){
  Theta <- Theta2Theta(Theta, k1=k, full1=full, check=check)
  A <- 1*(abs(Theta) > tol)
  graph <- igraph::graph_from_adjacency_matrix(
    A,
    mode = 'undirected',
    diag = FALSE
  )
  return(graph)
}
#' @rdname matrix2graph
#' @export
partialMatrixToGraph <- function(M){
  A <- !is.na(M)
  graph <- igraph::graph_from_adjacency_matrix(
    A,
    mode = 'undirected',
    diag = FALSE
  )
  return(graph)
}


#' Conversion between Huesler-Reiss parameter matrices
#' 
#' Converts between different matrices that parametrize the same
#' Huesler-Reiss distribution:
#' \eGamma, \eSigma, \eTheta, \eSigmaK, \eThetaK.
#' The \d1xd1 matrices \eSigmaK and \eThetaK can also be given/returned
#' as \dxd matrices with the kth row and column filled with zeros.
#' 
#' @inheritParams sharedParamsMatrixTransformations
#' 
#' @details
#' If `k`, `k1`, or `k2` is `NULL`, the corresponding `full*` argument is ignored.
#' 
#' @return
#' The desired parameter matrix corresponding to the specified inputs.
#' 
#' @family matrixTransformations
#' @rdname parameterMatrixConversion
#' @export
Gamma2Sigma <- function(Gamma, k=NULL, full=FALSE, check=TRUE){
  if(check){
    Gamma <- checkGamma(Gamma)
  }

  d <- nrow(Gamma)
  if(is.null(k)){
    ID <- diag(d) - matrix(1/d, d, d)
    Sigma <- ID %*% (-1/2 * Gamma) %*% ID
  } else{
    Sigma <- Sigma2Sigma(-1/2 * Gamma, k2=k, full2=full, check=FALSE)
  }

  if(check){
    Sigma <- ensure_matrix_symmetry_and_truncate_zeros(Sigma)
  }
  return(Sigma)
}

#' @rdname parameterMatrixConversion
#' @export
Gamma2Theta <- function(Gamma, k=NULL, full=FALSE, check=TRUE){
  Sigma <- Gamma2Sigma(Gamma, check=check)
  Theta <- Sigma2Theta(Sigma, k2=k, full2=full, check=FALSE)
  
  if(check){
    Theta <- ensure_matrix_symmetry_and_truncate_zeros(Theta)
  }
  
  return(Theta)
}

#' @rdname parameterMatrixConversion
#' @export
Sigma2Gamma <- function(Sigma, k=NULL, full=FALSE, check=TRUE){
  # The below transformation works for all dxd matrices that correspond to Gamma,
  # not just the one Sigma matrix -> only transform if not full
  if(!is.null(k) && !full){
    Sigma <- Sigma2Sigma(Sigma, k1=k, full1=full, check=check)
  } else if(check){
    Sigma <- checkSigma(Sigma, k, full)
  }

  d <- nrow(Sigma)
  oneVec <- makeOneVec(d)
  Gamma <- oneVec %*% t(diag(Sigma)) + diag(Sigma) %*% t(oneVec) - 2 * Sigma

  if(check){
    Gamma <- ensure_matrix_symmetry_and_truncate_zeros(Gamma)
  }
  return(Gamma)
}

#' @rdname parameterMatrixConversion
#' @export
Theta2Gamma <- function(Theta, k=NULL, full=FALSE, check=TRUE){
  Sigma <- Theta2Sigma(Theta, k1=k, full1=full, check=check)
  Gamma <- Sigma2Gamma(Sigma, check=FALSE)
  if(check){
    Gamma <- ensure_matrix_symmetry_and_truncate_zeros(Gamma)
  }
  return(Gamma)
}

#' @rdname parameterMatrixConversion
#' @export
Sigma2Theta <- function(Sigma, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
  # Convert input to Sigma (from Sigma^k) if necessary
  Sigma <- Sigma2Sigma(Sigma, k1=k1, full1=full1, check=check)
  
  Theta <- corpcor::pseudoinverse(Sigma)

  # Convert output to Theta^k if necessary
  if(!is.null(k2)){
    Theta <- Theta2Theta(Theta, k2=k2, full2=full2, check=FALSE)
  }

  if(check){
    Theta <- ensure_matrix_symmetry_and_truncate_zeros(Theta)
  }
  return(Theta)
}

#' @rdname parameterMatrixConversion
#' @export
Theta2Sigma <- function(Theta, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
  # Convert input to Theta (from Theta^k) if necessary
  Theta <- Theta2Theta(Theta, k1=k1, full1=full1, check=check)

  Sigma <- corpcor::pseudoinverse(Theta)
  
  # Convert output to Sigma^k if necessary
  if(!is.null(k2)){
    Sigma <- Sigma2Sigma(Sigma, k2=k2, full2=full2, check=FALSE)
  }
  
  if(check){
    Sigma <- ensure_matrix_symmetry_and_truncate_zeros(Sigma)
  }
  return(Sigma)
}

#' @rdname parameterMatrixConversion
#' @export
Theta2Theta <- function(Theta, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
  if(check){
    Theta <- checkTheta(Theta, k1, full1)
  }

  # Compute full Theta matrix
  if(full1 || is.null(k1)){
    ThetaFull <- Theta
  } else{
    d <- computeD(Theta, k1, full1)
    ThetaFull <- matrix(0, d, d)
    ThetaFull[-k1,-k1] <- Theta
  }
  if(!is.null(k1)){
    # Fill missing row/column s.t. rowSums/colSums are zero
    ThetaFull[,k1] <- -rowSums(ThetaFull)
    ThetaFull[k1,] <- ThetaFull[,k1]
    ThetaFull[k1,k1] <- -sum(ThetaFull[,k1])
  }
  
  # Compute return matrix
  if(check){
    ThetaFull <- ensure_matrix_symmetry_and_truncate_zeros(ThetaFull)
  }
  if(is.null(k2)){
    return(ThetaFull)
  }
  if(!full2){
    return(ThetaFull[-k2,-k2])
  }
  ThetaFull[k2,] <- 0
  ThetaFull[,k2] <- 0
  return(ThetaFull)
}

#' @rdname parameterMatrixConversion
#' @export
Sigma2Sigma <- function(Sigma, k1=NULL, k2=NULL, full1=FALSE, full2=FALSE, check=TRUE){
  if(check){
    Sigma <- checkSigma(Sigma, k1, full1)
  }

  d <- computeD(Sigma, k1, full1)

  # Compute full Sigma matrix
  if(full1 || is.null(k1)){
    SigmaFull <- Sigma
  } else{
    SigmaFull <- matrix(0, d, d)
    SigmaFull[-k1, -k1] <- Sigma
  }
  
  # Compute return matrix
  if(identical(k1, k2)){
    Sigma2 <- SigmaFull
    if(!is.null(k2) && !full2){
      Sigma2 <- Sigma2[-k2, -k2]
    }
  } else if(is.null(k2)){
    ID <- diag(d) - matrix(1/d, d, d)
    Sigma2 <- ID %*% SigmaFull %*% ID
  } else {
    Sigma2 <- (
      SigmaFull
      - matrix(SigmaFull[,k2], d, d, byrow = FALSE)
      - matrix(SigmaFull[,k2], d, d, byrow = TRUE)
      + SigmaFull[k2,k2]
    )

    if(!full2){
      Sigma2 <- Sigma2[-k2, -k2]
    }
  }

  # If requested, ensure numerical symmetry of return matrix
  if(check){
    Sigma2 <- ensure_matrix_symmetry_and_truncate_zeros(Sigma2)
  }
  return(Sigma2)
}

#' @details `Gamma2Gamma` only checks and returns the input.
#' 
#' @rdname parameterMatrixConversion
#' @export
Gamma2Gamma <- function(Gamma, check=TRUE){
  if(!check){
    return(Gamma) 
  }
  return(checkGamma(Gamma))
}

#' @param name1 Name of the input representation.
#' @param name2 Name of the output representation.
#' @param M Numeric matrix, \eGamma, \eSigma, or \eTheta.
#' @details `matrix2matrix` is a wrapper function that calls the corresponding
#' conversion function implied by `name1`, `name2`.
#' 
#' @rdname parameterMatrixConversion
#' @export
matrix2matrix <- function(
  M,
  name1=c('Gamma', 'Sigma', 'Theta')[1],
  name2=c('Gamma', 'Sigma', 'Theta')[1],
  k1=NULL,
  k2=NULL,
  full1=FALSE,
  full2=FALSE,
  check=TRUE
){
  ret <- if(name1 == 'Gamma'){
    if(name2 == 'Gamma'){
      Gamma2Gamma(M, check)
    } else if(name2 == 'Sigma'){
      Gamma2Sigma(M, k2, full2, check)
    } else if(name2 == 'Theta'){
      Gamma2Theta(M, k2, full2, check)
    }
  } else if(name1 == 'Sigma'){
    if(name2 == 'Gamma'){
      Sigma2Gamma(M, k1, full1, check)
    } else if(name2 == 'Sigma'){
      Sigma2Sigma(M, k1, k2, full1, full2, check)
    } else if(name2 == 'Theta'){
      Sigma2Theta(M, k1, k2, full1, full2, check)
    }
  } else if(name1 == 'Theta'){
    if(name2 == 'Gamma'){
      Theta2Gamma(M, k1, full1, check)
    } else if(name2 == 'Sigma'){
      Theta2Sigma(M, k1, k2, full1, full2, check)
    } else if(name2 == 'Theta'){
      Theta2Theta(M, k1, k2, full1, full2, check)
    }
  }
  if(is.null(ret)){
    stop('Invalid name1 or name2')
  }
  return(ret)
}


#' Create Gamma or Theta from vector
#'
#' Convert parameter vector `par` (upper triangular part of Gamma/Theta matrix)
#' to full Gamma/Theta, or vice versa.
#'
#' @param par Numeric vector with `d` elements.
#' Upper triangular part of a Gamma/Theta matrix.
#' @param allowMatrix If `TRUE` and `par` is already a matrix, return it as is.
#' @param allowNull If `TRUE` and `par` is NULL, return NULL.
#' @param zeroRowSums If `TRUE` the diagonal is set to (-1) times the rowSums.
#'
#' @return Numeric matrix \dxd. Full Gamma/Theta matrix corresponding to `par`.
#' 
#' @family matrixTransformations
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
  if(allowNull){
    return(NULL)
  }
  stop('Matrix `M` must not be NULL (unless allowNull=TRUE).')
  }
  if(allowVector && is.vector(M)){
  return(M)
  }
  return(upper.tri.val(M))
}


#' Transformation between \eChi and \eGamma
#'
#' Transforms between the extremal correlation \eChi and the variogram \eGamma.
#' Only valid for Huesler-Reiss distributions.
#' Done elementwise, no checks of the entire matrix structure are performed.
#'
#' @param chi Numeric vector or matrix with entries between 0 and 1.
#'
#' @return Numeric vector or matrix containing the implied \eGamma.
#'
#' @details
#' The formula for transformation from \eChi to \eGamma is element-wise 
#' \deqn{\Gamma = (2 \Phi^{-1}(1 - 0.5 \chi))^2,}
#' where \eqn{\Phi^{-1}} is the inverse of the standard normal distribution function.
#'
#' @family matrixTransformations
#' @rdname chi2Gamma
#' @export
chi2Gamma <- function(chi) {
  if (any(chi < 0) || any(chi > 1)) {
    stop("The argument chi must be between 0 and 1.")
  }
  Gamma <- (2 * stats::qnorm(1 - 0.5 * chi))^2
  return(Gamma)
}


#' @details
#' The formula for transformation from \eGamma to \eChi is element-wise 
#' \deqn{\chi = 2 - 2 \Phi(\sqrt{\Gamma} / 2),}{\chi = 2 - 2 \Phi(sqrt(\Gamma) / 2),}
#' where \eqn{\Phi} is the standard normal distribution function.
#'
#' @param Gamma Numeric vector or matrix with non-negative entries.
#' 
#' @return Numeric vector or matrix containing the implied \eChi.
#'
#' @rdname chi2Gamma
#' @export
Gamma2chi <- function(Gamma) {
  if(any(Gamma < 0)){
    stop("Entries of Gamma must be >= 0.")
  }
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
  d <- computeD(Gamma)

  if (d != 3) {
  stop("Gamma must be a 3 x 3 matrix.")
  }
  res <- 3 - V_HR(x = rep(1, times = 2), Gamma = Gamma[c(1, 2), c(1, 2)]) -
  V_HR(x = rep(1, times = 2), Gamma = Gamma[c(1, 3), c(1, 3)]) -
  V_HR(x = rep(1, times = 2), Gamma = Gamma[c(2, 3), c(2, 3)]) +
  V_HR(x = rep(1, times = 3), Gamma = Gamma)
  return(res)
}
