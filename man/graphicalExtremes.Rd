% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphicalExtremes-package.R
\docType{package}
\name{graphicalExtremes}
\alias{graphicalExtremes}
\alias{_PACKAGE}
\alias{graphicalExtremes-package}
\title{graphicalExtremes: Statistical methodology for graphical extreme value models.}
\description{
An implementation of the statistical methodology paper \localCiteT{eng2019} for sparse
multivariate extreme value models.
Includes exact simulation algorithms and statistical
inference methods for multivariate Pareto distributions on graphical structures.
Also contains implementations of statistical methods from
\localCiteT{eng2020}, \localCiteT{roe2021}, and \localCiteT{hen2022}.
}
\details{
The following global options are used by functions in the package.
Their values can be changed using \code{\link[base:options]{base::options()}}.
\describe{
\item{\code{"graphicalExtremes.mc.cores"}}{
The (maximal) number of cores to use in parallel tasks.
Will always be overwritten by 1 on Windows.
}
\item{\code{"graphicalExtremes.tol.small"}}{
The "small" tolerance is used in internal computations for values that should
mathematically be exactly \strong{equal to zero}, but deviate due to inherent
limitations of numerical computations. This value is used e.g. when checking
matrices for symmetry and definiteness.
In general, this value is used only as a "permissive" tolerance, in the sense
that if a value has to be positive, it is compared to actual zero, but if
it has to be zero, its absolute value is compared to this tolerance.
}
\item{\code{"graphicalExtremes.tol.large"}}{
The "large" tolerance is used for values that \strong{converge to zero}, but are
mathematically not supposed to be equal to zero. This value is used e.g.
when converting a precision matrix \eTheta to an adjacency matrix of a graph.
}
\item{\verb{"graphicalExtremes.default.alert}}{
The default alert function to be used in validity checks of Huesler-Reiss parameter matrix transformations.
Can be a function that takes an arbitrary number of strings as arguments (e.g. \code{cat()}, \code{stop()}),
\code{FALSE} to ignore the alerts, or \code{TRUE}/\code{NULL} to use the default function \code{warning()}.
}
}
}
\references{
\insertAllCited{}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/sebastian-engelke/graphicalExtremes}
  \item Report bugs at \url{https://github.com/sebastian-engelke/graphicalExtremes/issues}
}

}
\author{
\strong{Maintainer}: Sebastian Engelke \email{sebastian.engelke@unige.ch}

Authors:
\itemize{
  \item Adrien S. Hitz
  \item Nicola Gnecco
  \item Manuel Hentschel
}

}
\keyword{internal}
