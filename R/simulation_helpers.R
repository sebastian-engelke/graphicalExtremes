#' Simulate HR extremal functions
#'
#' Simulates the Huessler-Reiss extremal functions
#'
#' @inheritParams rmpareto
#' @param idx Integer. Index corresponding to the variable over which
#' the extremal function is simulated.
#' @param trend Numeric. Trend corresponding to the variable \code{idx}.
#' @param chol_mat Numeric matrix \eqn{d\times d}{d x d}.
#' Cholesky decomposition of the desired covariance matrix.
#' @return Numeric matrix \eqn{n\times d}{n x d}. Simulated data.
simu_px_HR <- function(n, idx, d, trend, chol_mat) {
  # check arguments
  if (length(idx)!=1) stop("Argument idx must be a scalar.")

  # function body
  d <- nrow(chol_mat)
  proc <- t(chol_mat)%*%matrix(rnorm(d*n), ncol=n) - trend
  proc <- exp(t(proc) - proc[idx,])
  return(proc)
}



#' Simulate logistic extremal functions
#'
#' Simulates logistic extremal functions
#'
#' @inheritParams simu_px_HR
#' @param idx Integer or numeric vector with \code{n} elements. Inde(x|ces) from
#' 1 to \code{d}, that determine which extremal function to simulate.
#' @param theta Numeric --- assume \eqn{0 < \theta < 1}.
#' @return Numeric matrix \eqn{n\times d}{n x d}. Simulated data.
simu_px_logistic <- function(n, idx, d, theta) {
  # check arguments
  if (length(idx) != 1 & length(idx) != n){
    stop("Argument idx must be a scalar or a vector with n entries")
  }



  # function body
  res <- matrix(1/gamma(1-theta)*(-log(runif(n*d)))^(-theta),
                      nrow=n, ncol=d)
  res[cbind(1:n,idx)] <- 1/gamma(1-theta)*rgamma(n,shape=1-theta)^(-theta)
  return(res/res[cbind(1:n,idx)])
}



#' Simulate negative logistic extremal functions
#'
#' Simulates negative logistic extremal functions
#'
#' @inheritParams simu_px_HR
#' @param idx Integer or numeric vector with \code{n} elements. Inde(x|ces) from
#' 1 to \code{d}, that determine which extremal function to simulate.
#' @param theta Numeric --- assume \eqn{\theta > 0}.
#' @return Numeric matrix \eqn{n\times d}{n x d}. Simulated data.
simu_px_neglogistic <- function(n, idx, d, theta) {
  # check arguments
  if (length(idx) != 1 & length(idx) != n){
    stop("Argument idx must be a scalar or a vector with n entries")
  }

  # function body
  res <- matrix(rweibull(n*d, shape=theta, scale=1/gamma(1+1/theta)),
                nrow=n, ncol=d)
  res[cbind(1:n,idx)] <- 1/gamma(1+1/theta)*rgamma(n,shape=1+1/theta)^(1/theta)
  return(res/res[cbind(1:n,idx)])
}



#' Simulate Dirichlet extremal functions
#'
#' Simulates Dirichlet extremal functions
#'
#' @inheritParams simu_px_HR
#' @param idx Integer or numeric vector with \code{n} elements. Inde(x|ces) from
#' 1 to \code{d}, that determine which extremal function to simulate.
#' @param theta Numeric vector of size \code{d}.
#' @return Numeric matrix \eqn{n\times d}{n x d}. Simulated data.
simu_px_dirichlet <- function(n, idx, d, alpha) {
  # check arguments
  if (length(idx) != 1 & length(idx) != n){
    stop("Argument idx must be a scalar or a vector with n entries")
  }

  # function body
  shape <- alpha
  shape[idx] <- alpha[idx] + 1
  shape.mat <- matrix(shape, nrow=d, ncol=n)
  rate.mat <- matrix(alpha, nrow=d, ncol=n)
  res <- t(matrix(rgamma(d*n, shape=shape.mat, rate=rate.mat), nrow=d, ncol=n))
  return(res/res[cbind(1:n,idx)])
}



#' Simulate HR extremal functions on a tree
#'
#' Simulates HR extremal functions on a tree
simu_px_tree_HR <- function(n, G.vec, A) {
  res <- exp(A %*% matrix(rnorm(length(G.vec)*n, mean= -G.vec/2, sd=sqrt(G.vec)), ncol=n))
  return(t(res))
}



### !!! Internal: simulates logistic extremal functions on a tree
simu_px_tree_logistic <- function(n, idx, nb.edges, theta, A) {
  stopifnot(length(idx)==1 || length(idx)==n)
  res       <- exp(A[[idx]] %*% log(matrix(1/gamma(1-theta)*(-log(runif(n*nb.edges)))^(-theta) /
                                             (1/gamma(1-theta)*rgamma(n*nb.edges,shape=1-theta)^(-theta)), ncol=n)))
  return(t(res))
}



### !!! Internal: simulates Dirichlet extremal functions on a tree
simu_px_tree_dirichlet <- function(n, alpha.start, alpha.end, A) {
  e = length(alpha.start)
  shape.start = matrix(alpha.start + 1, nrow=e, ncol=n)
  rate.start = matrix(alpha.start, nrow=e, ncol=n)
  shape.end = matrix(alpha.end, nrow=e, ncol=n)
  rate.end = matrix(alpha.end, nrow=e, ncol=n)
  sim.start = matrix(rgamma(e*n, shape=shape.start, rate=rate.start), nrow=e, ncol=n)
  sim.end = matrix(rgamma(e*n, shape=shape.end, rate=rate.end), nrow=e, ncol=n)
  res       <- exp(A %*% log(sim.end / sim.start))
  return(t(res))
}
