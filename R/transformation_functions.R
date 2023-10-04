

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
