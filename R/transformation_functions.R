
#' Data standardization to multivariate Pareto scale
#'
#' Transforms the `data` matrix empirically to the multivariate Pareto scale.
#'
#' @param data Numeric \nxd matrix, where `n` is the
#' number of observations and `d` is the dimension.
#' @param p Numeric between 0 and 1. Probability used for the quantile to
#' threshold the data.
#' @param na.rm Logical. If rows containing NAs should be removed.
#'
#' @return Numeric \eXTimesY{m}{d} matrix, where `m` is the number
#' of rows in the original `data` matrix that are above the threshold.
#'
#' @details
#' The columns of the `data` matrix are first transformed empirically to
#' standard Pareto distributions. Then, only the observations where at least
#' one component exceeds the `p`-quantile of the standard Pareto distribution
#' are kept. Those observations are finally divided by the `p`-quantile
#' of the standard Pareto distribution to standardize them to the multivariate Pareto scale.
#' 
#' If `na.rm` is `FALSE`, missing entries are left as such during the transformation of univariate marginals.
#' In the thresholding step, missing values are considered as `-Inf`.
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
#' 
#' @export
data2mpareto <- function(data, p, na.rm=FALSE) {
  # If specified, remove all rows that contain >=1 NA:
  if(na.rm){
    naInRow <- apply(is.na(data), 1, any)
    data <- data[!naInRow,,drop=FALSE]
  }
  if(nrow(data) == 0){
    return(data)
  }
  # Convert data (and quantile `p`) to std. Pareto marginals:
  dataPar <- matrix(
    1 / (1 - apply(data, 2, unif)), 
    nrow(data),
    ncol(data)
  )
  pPar <- 1 / (1 - p) 
  # Keep only rows with infty-norm > threshold:
  idx <- suppressWarnings(apply(dataPar, 1, max, na.rm=TRUE) > pPar)
  dataPar <- dataPar[idx,,drop=FALSE] / pPar
  return(dataPar)
}


#' Marginalize multivariate Pareto dataset
#'
#' Marginalize a multivariate Pareto dataset `data` with respect to the
#' variables in `set_indices`.
#'
#' @param data Numeric \nxd matrix. A dataset containing
#' observations following a multivariate Pareto distribution.
#' @param set_indices Numeric vector with at most `d` different elements in
#' 1, ..., `d`. The variables with respect to which to marginalize
#' the multivariate distribution.
#'
#' @return Numeric \eXTimesY{n}{m} matrix, where `m` is the length
#' of `set_indices`. Marginalized multivariate Pareto data.
#'
#' @keywords internal
mparetomargins <- function(data, set_indices) {
  data_sub <- data[, set_indices]
  idx <- which(apply(data_sub, 1, max) > 1)
  return(data[idx, set_indices])
}


#' Uniform margin
#'
#' Rescale the vector `x` empirically to uniform margin.
#'
#' @param x Numeric vector.
#' @param na.rm Logical. If TRUE, missing values are removed. If FALSE, missing values are kept as such.
#'
#' @return Numeric vector with entries rescaled to uniform margins
#'
#' @keywords internal
unif <- function(x, na.rm=FALSE) {
  if(na.rm){
    na.last <- NA
  } else{
    na.last <- 'keep'
  }
  rank(x, ties.method = "first", na.last=na.last) / (sum(!is.na(x)) + 1)
}
