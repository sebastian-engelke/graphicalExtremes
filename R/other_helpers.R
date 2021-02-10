#' Uniform margin
#'
#' Rescale the vector \code{x} empirically to uniform margin.
#'
#' @param x Numeric vector.
#'
#' @return Numeric vector with entries rescaled to uniform margins
#'
unif <- function(x) {
  rank(x) / (length(x) + 1)
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
dim_Gamma <- function(Gamma) {
  dimension <- dim(Gamma)

  if ((length(dimension) == 2) & (dimension[1] == dimension[2])) {
    dimension[1]
  } else {
    stop("Not a square matrix!")
  }
}



#' Select edges to add to a graph
#'
#' This function selects all possible edges that can be added to the \code{graph}
#' while still remaining in the class of block graphs.
#'
#' @inheritParams complete_Gamma
#'
#' @return Numeric vector.
#'
select_edges <- function(graph) {
  d <- igraph::vcount(graph)
  adj_mat <- igraph::as_adjacency_matrix(graph, sparse = FALSE) > 0


  sel.edges <- matrix(0, nrow = 0, ncol = 2)

  for (i in 1:(d - 1)) {
    for (j in (i + 1):d) {

      # set new_edge
      new_edge <- c(i, j)

      # check if new_edge is already in the graph
      is_already_edge <- adj_mat[i, j] | adj_mat[j, i]


      # if not, add it to the graph; else skip to the next
      if (!is_already_edge) {
        extended_graph <- igraph::add_edges(graph = graph, edges = new_edge)
      } else {
        next
      }

      # check if new graph is decomposable
      is_chordal <- igraph::is_chordal(extended_graph)$chordal

      # measure the length of the path from i to j in the old graph
      length_path <- length(as.vector(
        igraph::shortest_paths(graph, from = i, to = j)$vpath[[1]]
      ))

      if (is_chordal & length_path != 2) {
        sel.edges <- rbind(sel.edges, new_edge, deparse.level = 0)
      }
    }
  }

  return(sel.edges)
}



#' Set graphical parameters
#'
#' Set graphical parameters to \code{graph} which is an object from the
#' \code{igraph} package.
#'
#' @param graph Graph object from \code{igraph} package.
#'
#' @return Graph object from \code{igraph} package.
set_graph_parameters <- function(graph) {
  # set parameters
  igraph::V(graph)$color <- grDevices::adjustcolor(col = "#4477AA", alpha.f = 0.4)
  igraph::V(graph)$frame.color <- grDevices::adjustcolor(col = "#4477AA", alpha.f = 1)
  igraph::V(graph)$label.color <- "black"
  igraph::V(graph)$size <- 15
  igraph::E(graph)$width <- 2
  igraph::E(graph)$color <- "darkgrey"

  # return graph
  return(graph)
}


#' Censor dataset
#'
#' Censors each row of matrix \code{x} with vector \code{p}.
#'
#' @param x Numeric matrix \eqn{n \times d}{n x d}.
#' @param p Numeric vector with \eqn{d} elements.
#'
#' @return Numeric matrix \eqn{n \times d}{n x d}.
#'
censor <- function(x, p) {
  f2 <- function(x, p) {
    x_is_less <- x <= p
    y <- x
    y[x_is_less] <- p[x_is_less]
    return(y)
  }
  return(t(apply(x, 1, f2, p)))
}


is_eq <- function(a, b) {
  tol <- .Machine$double.eps^0.5
  abs(a - b) < tol
}

is_greater <- function(a, b) {
  tol <- .Machine$double.eps^0.5
  a - b > tol
}

is_less <- function(a, b) {
  is_greater(b, a)
}

is_leq <- function(a, b) {
  !(is_greater(a, b))
}

is_geq <- function(a, b) {
  !(is_less(a, b))
}

fast_diag <- function(y, M) {
  ## numeric_matrix numeric_matrix -> numeric_vector
  ## fast computation of diag(y %*% M %*% t(y))

  n <- nrow(y)
  sapply(1:n, function(i) {
    u <- y[i, , drop = FALSE]
    u %*% M %*% t(u)
  })
}
