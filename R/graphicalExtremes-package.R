#' graphicalExtremes: Statistical methodology for graphical extreme value models.
#'
#' An implementation of the statistical methodology paper \localCiteT{eng2019} for sparse
#' multivariate extreme value models.
#' Includes exact simulation algorithms and statistical 
#' inference methods for multivariate Pareto distributions on graphical structures.
#' Also contains implementations of statistical methods from
#' \localCiteT{eng2020}, \localCiteT{roe2021}, and \localCiteT{hen2022}.
#' 
#' The following global options are used by functions in the package.
#' Their values can be changed using [base::options()].
#' \describe{
#'  \item{`"graphicalExtremes.mc.cores"`}{
#'    The (maximal) number of cores to use in parallel tasks.
#'    Will always be overwritten by 1 on windows.
#'  }
#'  \item{`"graphicalExtremes.tol.small"`}{
#'    The "small" tolerance is used in internal computations for values that should
#'    mathematically be exactly **equal to zero**, but deviate due to inherent
#'    limitations of numerical computations. This value is used e.g. when checking
#'    matrices for symmetry and definiteness.
#'    In general, this value is used only as a "permissive" tolerance, in the sense
#'    that if a value has to be positive, it is compared to actual zero, but if
#'    it has to be zero, its absolute value is compared to this tolerance.
#'  }
#'  \item{`"graphicalExtremes.tol.large"`}{
#'    The "large" tolerance is used for values that **converge to zero**, but are
#'    mathematically not supposed to be equal to zero. This value is used e.g.
#'    when converting a precision matrix \eTheta to an adjacency matrix of a graph.
#'  }
#' }
#' 
#' 
#' @name graphicalExtremes
#
#' @importFrom Rdpack reprompt
#'
#' @references \insertAllCited{}
#' @aliases graphicalExtremes-package
#' @keywords internal
"_PACKAGE"
