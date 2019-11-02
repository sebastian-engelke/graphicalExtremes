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



#' Select edges to add to a graph
#'
#' This function selects all possible edges that can be added to graph
#' while still remaining in the class of block graphs.
#'
#' @inheritParams complete_Gamma
#'
#' @return Numeric vector.
#'
select_edges = function(graph){
  d = vcount(graph)
  sel.edges = matrix(0, nrow=0, ncol=2)
  for(i in 1:(d-1)) for(j in (i+1):d) if(is_chordal(add_edges(graph = graph, edges = c(i,j)))$chordal & length(as.vector(shortest_paths(graph, from=i, to=j)$vpath[[1]])) !=2) sel.edges = rbind(sel.edges,c(i,j))
  return(sel.edges)
}



#' Set graphical parameters
#'
#' Set graphical parameters to \code{graph} which is an object from the
#' \code{igraph} package.
#'
#' @inheritParams complete_Gamma
#'
#' @return Graph object from \code{igraph} package.
set_graph_parameters <- function(graph){
  # set parameters
  igraph::V(graph)$color <- adjustcolor(col = "#4477AA", alpha.f = 0.4)
  igraph::V(graph)$frame.color <- adjustcolor(col = "#4477AA", alpha.f = 1)
  igraph::V(graph)$label.color <- "black"
  igraph::V(graph)$size <- 15
  igraph::E(graph)$width <- 2
  igraph::E(graph)$color <- "darkgrey"

  # return graph
  return(graph)
}
