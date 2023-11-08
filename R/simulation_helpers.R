#' Simulate HR extremal functions
#'
#' Simulates the Huesler-Reiss extremal functions
#'
#' @inheritParams rmpareto
#' @param idx Integer. Index corresponding to the variable over which
#' the extremal function is simulated.
#' @param trend Numeric. Trend corresponding to the variable `idx`.
#' @param chol_mat Numeric matrix \dxd.
#' Cholesky decomposition of the desired covariance matrix.
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_HR <- function(n, d, idx, trend, chol_mat) {
  # check arguments
  if (length(idx) != 1) stop("Argument idx must be a scalar.")

  # function body
  d <- nrow(chol_mat)
  proc <- t(chol_mat) %*% matrix(stats::rnorm(d * n), ncol = n) - trend
  proc <- exp(t(proc) - proc[idx, ])
  return(proc)
}



#' Simulate logistic extremal functions
#'
#' Simulates logistic extremal functions
#'
#' @inheritParams simu_px_HR
#' @param idx Integer or numeric vector with `n` elements. Index/indices from
#' 1 to `d`, that determine which extremal function to simulate.
#' @param theta Numeric --- assume \eqn{0 < \theta < 1}.
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_logistic <- function(n, d, idx, theta) {
  # check arguments
  if (length(idx) != 1 & length(idx) != n) {
    stop("Argument idx must be a scalar or a vector with n entries")
  }

  # function body
  res <- matrix(1 / gamma(1 - theta) * (-log(stats::runif(n * d)))^(-theta),
    nrow = n, ncol = d
  )
  res[cbind(1:n, idx)] <- 1 / gamma(1 - theta) * stats::rgamma(n, shape = 1 - theta)^
    (-theta)
  return(res / res[cbind(1:n, idx)])
}



#' Simulate negative logistic extremal functions
#'
#' Simulates negative logistic extremal functions
#'
#' @inheritParams simu_px_HR
#' @param idx Integer or numeric vector with `n` elements. Index/indices from
#' 1 to `d`, that determine which extremal function to simulate.
#' @param theta Numeric --- assume \eqn{\theta > 0}.
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_neglogistic <- function(n, d, idx, theta) {
  # check arguments
  if (length(idx) != 1 & length(idx) != n) {
    stop("Argument idx must be a scalar or a vector with n entries")
  }

  # function body
  res <- matrix(stats::rweibull(n * d, shape = theta, scale = 1 /
    gamma(1 + 1 / theta)), nrow = n, ncol = d)
  res[cbind(1:n, idx)] <- 1 / gamma(1 + 1 / theta) *
    stats::rgamma(n, shape = 1 + 1 / theta)^(1 / theta)
  return(res / res[cbind(1:n, idx)])
}



#' Simulate Dirichlet extremal functions
#'
#' Simulates Dirichlet extremal functions
#'
#' @inheritParams simu_px_HR
#' @param idx Integer or numeric vector with `n` elements. Index/indices from
#' 1 to `d`, that determine which extremal function to simulate.
#' @param alpha Numeric vector of size `d`.
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_dirichlet <- function(n, d, idx, alpha) {
  # check arguments
  if (length(idx) != 1 & length(idx) != n) {
    stop("Argument idx must be a scalar or a vector with n entries")
  }

  # function body
  shape <- alpha
  shape[idx] <- alpha[idx] + 1
  shape.mat <- matrix(shape, nrow = d, ncol = n)
  rate.mat <- matrix(alpha, nrow = d, ncol = n)
  res <- t(matrix(stats::rgamma(d * n, shape = shape.mat, rate = rate.mat),
    nrow = d, ncol = n
  ))
  return(res / res[cbind(1:n, idx)])
}



#' Simulate HR extremal functions on a tree
#'
#' Simulates the Huesler-Reiss extremal functions on a tree
#'
#' @inheritParams rmpareto
#' @param Gamma_vec Numeric vector with `d-1` elements, where `d` is the
#' number of nodes in the tree (and `d-1` is the number of edges).
#' @param A Numeric \dxd1 matrix; the rows represent the nodes
#' in the tree, the columns represent the edges. For a fixed node
#' \eqn{k = 1, \dots, d}{k = 1, ..., d}, each entry \eqn{(i, j)} is
#' equal to 1 if the edge in position `j` is on the directed path from node
#' `k` to node `i` in the polytree rooted at node `k`.
#'
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_tree_HR <- function(n, Gamma_vec, A) {
  res <- exp(A %*%
    matrix(stats::rnorm(length(Gamma_vec) * n,
      mean = -Gamma_vec / 2,
      sd = sqrt(Gamma_vec)
    ),
    ncol = n
    ))
  return(t(res))
}



#' Simulate logistic extremal functions on a tree
#'
#' Simulates logistic extremal functions on a tree
#'
#' @inheritParams rmpareto
#' @param theta Numeric vector with 1 or `d-1` elements.
#' Assume that the entry are such that \eqn{0 < \theta < 1}.
#' @param A Numeric \dxd1 matrix; the rows represent the
#' nodes in the tree, the columns represent the edges. For a fixed node
#' \eqn{k = 1, \dots, d}{k = 1, ..., d}, each entry \eqn{(i, j)} is
#' equal to 1 if the edge in position `j` is on the directed path from node
#' `k` to node `i` in the polytree rooted at node `k`.
#'
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_tree_logistic <- function(n, theta, A) {
  # define number of edges
  d <- nrow(A)
  nb.edges <- d - 1

  # check arguments
  if (length(theta) != 1 & length(theta) != d - 1) {
    stop("Argument theta must be a vector with 1 or d - 1 entries")
  }




  # function body
  res <- exp(A %*%
    log(matrix(1 / gamma(1 - theta) *
      (-log(stats::runif(n * nb.edges)))^(-theta) /
      (1 / gamma(1 - theta) *
        stats::rgamma(n * nb.edges, shape = 1 - theta)^
          (-theta)), ncol = n)))
  return(t(res))
}



#' Simulate Dirichlet extremal functions on a tree
#'
#' Simulates Dirichlet extremal functions on a tree
#'
#' @inheritParams rmpareto
#' @param alpha.start Numeric vector with `d-1` elements, where `d` is
#' the number of nodes in the tree (and `d-1` is the number of edges).
#' @param alpha.end Numeric vector with `d-1` elements, where `d` is
#' the number of nodes in the tree (and `d-1` is the number of edges).
#' @param A Numeric \dxd1 matrix; the rows represent the
#' nodes in the tree, the columns represent the edges. For a fixed node
#' \eqn{k = 1, \dots, d}{k = 1, ..., d}, each entry \eqn{(i, j)} is
#' equal to 1 if the edge in position `j` is on the directed path from node
#' `k` to node `i` in the polytree rooted at node `k`.
#'
#' @return Numeric \nxd matrix. Simulated data.
#'
#' @keywords internal
simu_px_tree_dirichlet <- function(n, alpha.start, alpha.end, A) {
  e <- length(alpha.start)
  shape.start <- matrix(alpha.start + 1, nrow = e, ncol = n)
  rate.start <- matrix(alpha.start, nrow = e, ncol = n)
  shape.end <- matrix(alpha.end, nrow = e, ncol = n)
  rate.end <- matrix(alpha.end, nrow = e, ncol = n)
  sim.start <- matrix(stats::rgamma(e * n, shape = shape.start, rate = rate.start),
    nrow = e, ncol = n
  )
  sim.end <- matrix(stats::rgamma(e * n, shape = shape.end, rate = rate.end),
    nrow = e, ncol = n
  )
  res <- exp(A %*% log(sim.end / sim.start))
  return(t(res))
}
