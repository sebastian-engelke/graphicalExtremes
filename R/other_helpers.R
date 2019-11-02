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



#' Is Gamma square matrix?
#'
#' Check if Gamma matrix is square matrix. If so, return the dimension. Else,
#' raise an error.
#'
#' @param Gamma Numeric matrix. Matrix representing the variogram of an HR
#' distribution.
#'
#' @return Numeric. The dimension of the matrix (number of rows and columns, if
#' the matrix is symmetric). Else, raises an error.
#'
dim_Gamma <- function(Gamma){
  dimension <- dim(Gamma)

  if ((length(dimension) == 2) & (dimension[1] == dimension[2])){
    dimension[1]
  } else {
    stop("Not a square matrix!")
  }
}
