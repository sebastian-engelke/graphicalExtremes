% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate_examples.R
\name{generate_random_model}
\alias{generate_random_model}
\title{Generate random Hüsler--Reiss Models}
\usage{
generate_random_model(d, graph_type = "tree", ...)
}
\arguments{
\item{d}{Number of vertices in the graph}

\item{graph_type}{\code{"tree"}, \code{"block"}, \code{"decomposable"}, or \code{"general"}}

\item{...}{Further arguments passed to functions generating the graph and Gamma matrix}
}
\description{
Generates a random connected graph and Gamma matrix with conditional independence
structure corresponding to that graph.
}
\examples{
set.seed(1)
d <- 12

generate_random_model(d, 'tree')
generate_random_model(d, 'block')
generate_random_model(d, 'decomposable')
generate_random_model(d, 'general')

}
\seealso{
Other Example generations: 
\code{\link{generate_random_chordal_graph}()},
\code{\link{generate_random_connected_graph}()},
\code{\link{generate_random_graphical_Gamma}()},
\code{\link{generate_random_integer_Gamma}()},
\code{\link{generate_random_spd_matrix}()},
\code{\link{generate_random_tree}()}
}
\concept{Example generations}
