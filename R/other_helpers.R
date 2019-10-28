#' Uniform margin
#'
#' Rescale the vector \code{x} to uniform margin.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector with entries rescaled to uniform margins
#'
unif <- function(x){
  rank(x)/(length(x)+1)
}
