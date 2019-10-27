#' chi.est
#'
#' Estimates the chi coefficient empirically.
#'
#' The \eqn{\chi} coefficient is a scalar coefficient that represents
#' the extremal correlation between two variables.
#'
#' @param data Numeric matrix \eqn{n \times 2}{n x 2}. A data matrix
#' containing two variables.
#' @param u Numeric, between 0 and 1. It is the probability threshold used to
#' compute the \eqn{\chi} coefficient.
#' @param pot Boolean. if TRUE, then pot-type estimation of EC is used.
#' By default, \code{pot = FALSE}.
#'
#' @return Numeric. The empirical \eqn{\chi} coefficient between the 2 variables
#' in \code{data}.
#'
chi.est <- function(data, u, pot=FALSE){
  if (NCOL(data) > 2){
    warning("The data matrix should contain only two columns.")
  }

  data <- na.omit(data)
  n <- nrow(data)
  data <-  apply(data, 2, unif)
  rowmax <- apply(data, 1, max)
  rowmin <- apply(data, 1, min)
  eps <- .Machine$double.eps^0.5
  qlim2 <- c(min(rowmax) + eps, max(rowmin) - eps)

  qlim <- qlim2
  nq <- length(u)
  cu <- cbaru <- numeric(nq)
  for (i in 1:nq) cu[i] <- mean(rowmax < u[i])
  for (i in 1:nq) cbaru[i] <- mean(rowmin > u[i])
  if(pot) chiu <- cbaru / (1-u)
  if(!pot) chiu <- 2 - log(cu)/log(u)
  return(chiu)
}



#' chi_mat
#'
#' Estimates empirically the extremal \eqn{\chi} coefficient
#' given the dataset \code{data}.
#'
#' @param data Numeric matrix \eqn{n \times d}{n x d}. A data matrix
#' containing \eqn{d} variables.
#' @inheritParams chi.est
#'
#' @return Numeric matrix \eqn{d\times d}{d x d}. The matrix containing the
#' bivariate extremal coefficientes \eqn{\chi_{ij}}, for \eqn{i, j = 1, ..., d}.
#'
chi_mat <- function(data, u, pot=FALSE){
  d <- ncol(data)
  res <- as.matrix(expand.grid(1:d,1:d))
  res <- res[res[,1]>res[,2],,drop=FALSE]
  chi <- apply(res, 1, function(x){
    chi.est(cbind(data[,x[1]], data[,x[2]]), u=u, pot=pot)
  })[1]
  chi.mat <- matrix(NA, ncol=d, nrow=d)
  chi.mat[res] <- chi
  chi.mat[res[,2:1]] <- chi
  diag(chi.mat) <- 1

  return(chi.mat)
}
