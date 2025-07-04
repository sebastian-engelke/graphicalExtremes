% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimation_structure.R
\name{eglearn}
\alias{eglearn}
\title{Learning extremal graph structure}
\usage{
eglearn(
  data,
  p = NULL,
  rholist = c(0.1, 0.15, 0.19, 0.205),
  reg_method = c("ns", "glasso"),
  complete_Gamma = FALSE
)
}
\arguments{
\item{data}{Numeric \nxd matrix, where \code{n} is the
number of observations and \code{d} is the dimension.}

\item{p}{Numeric between 0 and 1 or \code{NULL}. If \code{NULL} (default),
it is assumed that the \code{data} are already on multivariate Pareto scale. Else,
\code{p} is used as the probability in the function \code{\link[=data2mpareto]{data2mpareto()}}
to standardize the \code{data}.}

\item{rholist}{Numeric vector of non-negative regularization parameters
for the lasso.
Default is \code{rholist = c(0.1, 0.15, 0.19, 0.205)}.
For details see \code{glasso::glassopath()}.}

\item{reg_method}{One of \verb{"ns", "glasso"}, for neighborhood selection and
graphical lasso, respectively.
Default is \code{reg_method = "ns"}.
For details see \insertCite{meins2006;textual}{graphicalExtremes},
\insertCite{friedman2008;textual}{graphicalExtremes}.}

\item{complete_Gamma}{Whether you want to try fto complete Gamma matrix.
Default is \code{complete_Gamma = FALSE}.}
}
\value{
List made of:
\item{\code{graph}}{
A list of \code{\link[igraph:make_graph]{igraph::graph}} objects representing the
fitted graphs for each \code{rho} in \code{rholist}.
}
\item{\code{Gamma}}{
A list of numeric estimated \dxd
variogram matrices \eGamma corresponding to the fitted graphs,
for each \code{rho} in \code{rholist}. If \code{complete_Gamma = FALSE} or the
underlying graph is not connected, it returns \code{NULL}.
}
\item{\code{rholist}}{
The list of penalty coefficients.
}
\item{\code{graph_ic}}{
A list of \code{\link[igraph:make_graph]{igraph::graph}} objects
representing the optimal graph
according to the \code{aic}, \code{bic}, and \code{mbic} information criteria.
If \code{reg_method = "glasso"}, it returns a list of \code{NULL}.
}
\item{\code{Gamma_ic}}{
A list of numeric \dxd estimated
variogram matrices \eGamma corresponding
to the \code{aic}, \code{bic}, and \code{mbic} information criteria.
If \code{reg_method = "glasso"}, \code{complete_Gamma = FALSE}, or the underlying
graph is not connected, it returns a list of \code{NULL}.
}
}
\description{
Following the methodology from \insertCite{eng2022a;textual}{graphicalExtremes},
fits an extremal graph structure using the neighborhood selection approach
(see \insertCite{meins2006;textual}{graphicalExtremes}) or graphical lasso
(see \insertCite{friedman2008;textual}{graphicalExtremes}).
}
\examples{
set.seed(2)
m <- generate_random_model(d=6)
y <- rmpareto(n=500, par=m$Gamma)
ret <- eglearn(y)

}
\references{
\insertAllCited{}
}
\seealso{
Other structure estimation methods:
\code{\link{data2mpareto}()},
\code{\link{eglatent}()},
\code{\link{emst}()},
\code{\link{fit_graph_to_Theta}()}
}
\concept{structureEstimation}
