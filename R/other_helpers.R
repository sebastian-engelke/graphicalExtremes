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
#' @keywords internal
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
#' This function selects all possible edges that can be added to the `graph`
#' while still remaining in the class of block graphs.
#'
#' @inheritParams complete_Gamma
#'
#' @return Numeric vector.
#'
#' @keywords internal
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


graphs_equal <- function(g1, g2) {
  ## graph graph -> boolean
  ## produce true if two graphs have same edges
  
  # Return early if graph sizes are different
  if(igraph::vcount(g1) != igraph::vcount(g2)){
    return(FALSE)
  }
  # Compare adjacency matrices
  A1 <- igraph::as_adjacency_matrix(g1, sparse = FALSE)
  A2 <- igraph::as_adjacency_matrix(g2, sparse = FALSE)
  return(all(A1 == A2))
}
