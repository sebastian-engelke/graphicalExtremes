
#' Transformation of \eqn{\Gamma} matrix to graph object
#'
#' Transforms \code{Gamma} matrix to an \code{igraph} object for
#' the corresponding Huesler--Reiss extremal graphical model,
#' and plots it (optionally).
#'
#' @param Gamma Numeric \eqn{d \times d}{d x d} variogram matrix.
#' @param to_plot Logical. If \code{TRUE} (default), it plots the resulting
#' graph.
#' @param ... Graphical parameters for the \code{\link[igraph]{plot.igraph}}
#' function of the package \code{igraph}.
#'
#' @return Graph object from \code{igraph} package. An undirected graph.
#'
#' @details
#' The variogram uniquely determines the extremal graph structure of the
#' corresponding Huesler--Reiss distribution. The conditional independencies
#' can be identified from the inverses of the matrices \eqn{\Sigma^{(k)}}{\Sigma^(k)}
#' defined in equation (10) in \insertCite{eng2019;textual}{graphicalExtremes}.
#'
#' @examples
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' Gamma2graph(Gamma, to_plot = TRUE)
#' @references
#'  \insertAllCited{}
#'
#' @export
Gamma2graph <- function(Gamma, to_plot = TRUE, ...) {
  null.mat <- matrix(0, nrow = nrow(Gamma), ncol = ncol(Gamma))
  for (i in 1:nrow(Gamma)) {
    null.mat[-i, -i] <- null.mat[-i, -i] +
      (abs(solve(Gamma2Sigma(Gamma, i))) < 1e-6)
  }
  graph <- igraph::graph_from_adjacency_matrix(null.mat == 0,
    diag = FALSE,
    mode = "undirected"
  )

  if (to_plot) {
    igraph::plot.igraph(graph, ...)
  }
  return(graph)
}



#' Data standardization to multivariate Pareto scale
#'
#' Transforms the \code{data} matrix empirically to the multivariate Pareto scale.
#'
#' @param data Numeric matrix of size \eqn{n\times d}{n x d}, where \eqn{n} is the
#' number of observations and \eqn{d} is the dimension.
#' @param p Numeric between 0 and 1. Probability used for the quantile to
#' threshold the data.
#'
#' @return Numeric matrix \eqn{m \times d}{m x d}, where \eqn{m} is the number
#' of rows in the original \code{data} matrix that are above the threshold.
#'
#' @details
#' The columns of the \code{data} matrix are first transformed empirically to
#' standard Pareto distributions. Then, only the observations where at least
#' one component exceeds the \code{p}-quantile of the standard Pareto distribution
#' are kept. Those observations are finally divided by the \code{p}-quantile
#' of the standard Pareto distribution to standardize them to the multivariate Pareto scale.
#'
#' @examples
#' n <- 20
#' d <- 4
#' p <- .8
#' G <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#'
#' set.seed(123)
#' my_data <- rmstable(n, "HR", d = d, par = G)
#' data2mpareto(my_data, p)
#' @export
data2mpareto <- function(data, p) {
  xx <- 1 / (1 - apply(data, 2, unif))
  q <- stats::quantile(xx, p)
  idx <- which(apply(xx, 1, max) > q)
  return(xx[idx, ] / q)
}


#' Transformation of \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix to \eqn{\Gamma} matrix
#'
#' Transforms the \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix from the definition of a
#' Huesler--Reiss distribution to the corresponding \eqn{\Gamma} matrix.
#'
#' @param S Numeric \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)} covariance matrix \eqn{\Sigma^{(k)}}{\Sigma^(k)}
#' from the definition of a Huesler--Reiss distribution.
#' Numeric \eqn{d \times d}{d x d} covariance matrix if \code{full = TRUE}, see \code{full}
#' parameter.
#' @param k Integer between \code{1} (the default value) and \code{d}.
#' Indicates which matrix \eqn{\Sigma^{(k)}}{\Sigma^(k)} is represented by \code{S}.
#' @param full Logical. If true, then the \code{k}th row and column in \eqn{\Sigma^{(k)}}{\Sigma^(k)}
#' are included and the function returns a \eqn{d \times d}{d x d} matrix.
#' By default, \code{full = FALSE}.
#'
#' @details
#' For any \code{k} from \code{1} to \code{d},
#' the \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix of size \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)}
#' in the definition of a
#' Huesler--Reiss distribution can be transformed into a the
#' corresponding \eqn{d \times d}{d x d} \eqn{\Gamma} matrix.
#' If \code{full = TRUE}, then \eqn{\Sigma^{(k)}}{\Sigma^(k)} must be a \eqn{d \times d}{d x d}
#' matrix with \code{k}th row and column
#' containing zeros. For details see \insertCite{eng2019;textual}{graphicalExtremes}.
#' This is the inverse of function of \code{\link{Gamma2Sigma}}.
#'
#' @return Numeric \eqn{d \times d}{d x d} \eqn{\Gamma} matrix.
#'
#' @examples
#' Sigma1 <- rbind(
#'   c(1.5, 0.5, 1),
#'   c(0.5, 1.5, 1),
#'   c(1, 1, 2)
#' )
#' Sigma2Gamma(Sigma1, k = 1, full = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
Sigma2Gamma <- function(S, k = 1, full = FALSE) {
  # complete S
  if (!full) {
    d <- NROW(S)
    S_full <- rbind(rep(0, d + 1), cbind(rep(0, d), S))

    shuffle <- 1:(d + 1)
    shuffle[shuffle <= k] <- shuffle[shuffle <= k] - 1
    shuffle[1] <- k
    shuffle <- order(shuffle)

    S_full <- S_full[shuffle, shuffle]
  } else {
    S_full <- S
  }

  # compute Gamma
  One <- rep(1, times = ncol(S_full))
  D <- diag(S_full)
  Gamma <- One %*% t(D) + D %*% t(One) - 2 * S_full

  return(Gamma)
}



#' Transformation of \eqn{\Gamma} matrix to \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix
#'
#' Transforms the \code{Gamma} matrix from the definition of a
#' Huesler--Reiss distribution to the corresponding \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix.
#'
#'
#' @param Gamma Numeric \eqn{d \times d}{d x d} variogram matrix.
#' @param k Integer between \code{1} (the default value) and \code{d}.
#' Indicates which matrix \eqn{\Sigma^{(k)}}{\Sigma^(k)} should be produced.
#' @param full Logical. If true, then the \code{k}th row and column in \eqn{\Sigma^{(k)}}{\Sigma^(k)}
#' are included and the function returns a \eqn{d \times d}{d x d} matrix.
#' By default, \code{full = FALSE}.
#'
#' @details
#' Every \eqn{d \times d}{d x d} \code{Gamma} matrix in the definition of a
#' Huesler--Reiss distribution can be transformed into a
#' \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)} \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix,
#' for any \code{k} from \code{1} to \code{d}. The inverse of \eqn{\Sigma^{(k)}}{\Sigma^(k)}
#' contains the graph structure corresponding to the Huesler--Reiss distribution
#' with parameter matrix \code{Gamma}. If \code{full = TRUE}, then \eqn{\Sigma^{(k)}}{\Sigma^(k)}
#' is returned as a \eqn{d \times d}{d x d} matrix with additional \code{k}th row and column
#' that contain zeros. For details see \insertCite{eng2019;textual}{graphicalExtremes}.
#' This is the inverse of function of \code{\link{Sigma2Gamma}}.
#'
#' @return Numeric \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix of size \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)} if
#' \code{full = FALSE}, and of size \eqn{d \times d}{d x d} if \code{full = TRUE}.
#'
#' @examples
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#' Gamma2Sigma(Gamma, k = 1, full = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
Gamma2Sigma <- function(Gamma, k = 1, full = FALSE) {
  d <- ncol(Gamma)
  if (full) {
    1 / 2 * (matrix(rep(Gamma[, k], d), ncol = d, nrow = d) +
      t(matrix(rep(Gamma[, k], d), ncol = d, nrow = d)) - Gamma)
  } else {
    1 / 2 * (matrix(rep(Gamma[-k, k], d - 1), ncol = d - 1, nrow = d - 1) +
      t(matrix(rep(Gamma[-k, k], d - 1), ncol = d - 1, nrow = d - 1)) - Gamma[-k, -k])
  }
}


#' Transformation of \eqn{\Gamma} matrix to \eqn{\Theta} matrix
#'
#' Transforms the \code{Gamma} matrix from the definition of a
#' Huesler--Reiss distribution to the corresponding \eqn{\Theta} matrix.
#'
#'
#' @param Gamma Numeric \eqn{d \times d}{d x d} variogram matrix.
#'
#' @details
#' Every \eqn{d \times d}{d x d} \code{Gamma} matrix in the definition of a
#' Huesler--Reiss distribution can be transformed into a
#' \eqn{d \times d}{d x d} \eqn{\Theta} matrix, which
#' contains the graph structure corresponding to the Huesler--Reiss distribution
#' with parameter matrix \code{Gamma}.
#'
#' @return Numeric \eqn{\Sigma^{(k)}}{\Sigma^(k)} matrix of size \eqn{(d - 1) \times (d - 1)}{(d - 1) x (d - 1)} if
#' \code{full = FALSE}, and of size \eqn{d \times d}{d x d} if \code{full = TRUE}.
#'
#' @examples
#' Gamma <- cbind(
#'   c(0, 1.5, 1.5, 2),
#'   c(1.5, 0, 2, 1.5),
#'   c(1.5, 2, 0, 1.5),
#'   c(2, 1.5, 1.5, 0)
#' )
#' Gamma2Sigma(Gamma, k = 1, full = FALSE)
#' @references
#'  \insertAllCited{}
#'
#' @export
Gamma2Theta <- function(Gamma, k=NULL) {
  if(is.null(k)){
    d <- ncol(Gamma)
    ID <- diag(d) - matrix(1/d, d, d)
    S <- ID %*% (-1/2 * Gamma) %*% ID
    Theta <- corpcor::pseudoinverse(S)
    return(Theta)
  } else{
    Sigma_k <- Gamma2Sigma(Gamma, k)
    Theta_k <- chol2inv(chol(Sigma_k))
    return(Theta_k)
  }
}

#' Transformation of \eqn{\Gamma} matrix to \eqn{\Theta} matrix
#'
#' @export
Theta2Gamma <- function(Theta, k=NULL) {
  if(is.null(k)){
    d <- nrow(Theta)
    id <- matrix(1, d)
    S <- corpcor::pseudoinverse(Theta)
    dS <- diag(S)
    return(id %*% t(dS) + dS %*% t(id) - 2*S)
  } else{
    Sigma_k <- chol2inv(chol(Theta))
    return(Sigma2Gamma(Sigma_k, k))
  }
}


#' Create \eqn{\Gamma} from vector
#'
#' This function takes the parameters in the vector \code{par}
#' (upper triangular Gamma matrix) and returns full Gamma.
#'
#' @param par Numeric vector with \eqn{d} elements.
#' Upper triangular part of a Gamma matrix.
#'
#' @return Numeric matrix \eqn{d \times d}{d x d}. Full Gamma matrix.
#'
#' @keywords internal
par2Gamma <- function(par) {
  d <- 1 / 2 + sqrt(1 / 4 + 2 * length(par))
  if (round(d) != d) {
    stop("The length of par does not agree with any square matrix.")
  }
  G <- matrix(0, nrow = d, ncol = d)
  G[upper.tri(G)] <- par
  return(G + t(G))
}



#' Extract upper triangular part of \eqn{\Gamma}
#'
#' This function returns a vector containing the upper triangular part
#' of the matrix \code{Gamma}. If \code{Gamma} is already a vector, it returns
#' it as it is.
#'
#' @param Gamma Numeric \eqn{d \times d}{d x d} variogram matrix.
#'
#' @return Numeric vector with \eqn{d} elements.
#' The upper triangular part of the given \code{Gamma} matrix.
#'
#' @keywords internal
Gamma2par <- function(Gamma) {
  if (is.matrix(Gamma)) {
    return(Gamma[upper.tri(Gamma)])
  } else {
    return(Gamma)
  }
}



#' Transformation of extremal correlation \eqn{\chi} to the Huesler--Reiss variogram \eqn{\Gamma}
#'
#' Transforms the extremal correlation \eqn{\chi} into the \code{Gamma} matrix
#' from the definition of a Huesler--Reiss
#' distribution.
#'
#' @param chi Numeric or matrix, with entries
#' between 0 and 1.
#'
#' @return Numeric or matrix. The \eqn{\Gamma} parameters in the Huesler--Reiss
#' distribution.
#'
#' @details
#' The formula for transformation from \code{chi} to \eqn{\Gamma} that is applied element-wise is
#' \deqn{\Gamma = (2 \Phi^{-1}(1 - 0.5 \chi))^2,}
#' where \eqn{\Phi^{-1}} is the inverse of the standard normal distribution function.
#' This is the inverse of \code{\link{Gamma2chi}}.
#'
#' @export
chi2Gamma <- function(chi) {
  if (any(is_less(chi, 0) | is_greater(chi, 1))) {
    stop("The argument chi must be between 0 and 1.")
  }
  Gamma <- (2 * stats::qnorm(1 - 0.5 * chi))^2
  return(Gamma)
}


#' Transformation of the Huesler--Reiss variogram \eqn{\Gamma} to extremal correlation \eqn{\chi}
#'
#' Transforms the \code{Gamma} matrix from the definition of a Huesler--Reiss
#' distribution into the corresponding extremal correlation \eqn{\chi}.
#'
#' @details
#' The formula for transformation from \code{Gamma} to \eqn{\chi} that is applied element-wise is
#' \deqn{\chi = 2 - 2 \Phi(\sqrt{\Gamma} / 2),}{\chi = 2 - 2 \Phi(sqrt(\Gamma) / 2),}
#' where \eqn{\Phi} is the standard normal distribution function.
#' This is the inverse of \code{\link{chi2Gamma}}.
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




#' Compute theoretical \eqn{\chi} in 3D
#'
#' Computes the theoretical \eqn{\chi} coefficient in 3 dimensions.
#'
#' @param Gamma Numeric matrix \eqn{3\times 3}{3 x 3}.
#'
#' @return The 3-dimensional \eqn{\chi} coefficient, i.e.,
#' the extremal correlation coefficient for the HR distribution. Note that
#' \eqn{0 \leq \chi \leq 1}.
#'
#' @keywords internal
Gamma2chi_3D <- function(Gamma) {
  d <- dim_Gamma(Gamma)

  if (d != 3) {
    stop("Gamma must be a 3 x 3 matrix.")
  }
  res <- 3 - V_HR(x = rep(1, times = 2), par = Gamma2par(Gamma[c(1, 2), c(1, 2)])) -
    V_HR(x = rep(1, times = 2), par = Gamma2par(Gamma[c(1, 3), c(1, 3)])) -
    V_HR(x = rep(1, times = 2), par = Gamma2par(Gamma[c(2, 3), c(2, 3)])) +
    V_HR(x = rep(1, times = 3), par = Gamma2par(Gamma))
  return(res)
}



#' Marginalize multivariate Pareto dataset
#'
#' Marginalize a multivariate Pareto dataset \code{data} with respect to the
#' variables in \code{set_indices}.
#'
#' @param data Numeric matrix \eqn{n\times d}{n x d}. A dataset containing
#' observations following a multivariate Pareto distribution.
#' @param set_indices Numeric vector with at most \eqn{d} different elements in
#' 1, ..., \eqn{d}. The variables with respect to which to marginalize
#' the multivariate distribution.
#'
#' @return Numeric matrix \eqn{n\times m}{n x m}, where \eqn{m} is the length
#' of \code{set_indices}. Marginalized multivariate Pareto data.
#'
#' @keywords internal
mparetomargins <- function(data, set_indices) {
  data_sub <- data[, set_indices]
  idx <- which(apply(data_sub, 1, max) > 1)
  return(data[idx, set_indices])
}
