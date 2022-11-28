#' Check input graph
#'
#' Checks that the input graph is a valid graph for an extremal graphical model.
#' If necessary, converts the graph to an undirected graph.
#' Removes vertex labels if present.
#'
#' @param graph An [igraph::graph] object.
#' @param graph_type `"general"`, `"decomposable"`, `"block"`, `"tree"`. The required type of graph.
#' @param check_connected Whether to check if the graph is connected.
#' @param nVertcies The number of vertices required in the graph.
#'
#' @return The given `graph`, if necessary converted to undirected.
#' If the graph is not valid an error is thrown.
#'
#' @family Input checks
check_graph <- function(
  graph,
  graph_type = c('general', 'decomposable', 'block', 'tree'),
  check_connected = TRUE,
  nVertices = NULL
){
  graph_type <- match.arg(graph_type)

  # check that graph is actually a graph object
  if(!igraph::is.igraph(graph)){
    stop('The given object is not an igraph-graph.')
  }

  # set up main variables
  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # check graph size
  if(!is.null(nVertices) && d != nVertices){
    stop('The given graph does not have the correct number of vertices.')
  }

  # check if it is directed
  if (igraph::is_directed(graph)) {
    graph <- igraph::as.undirected(graph)
  }

  if(!is.null(igraph::vertex_attr(graph)[['name']])){
    graph <- igraph::remove.vertex.attribute(graph, 'name')
  }

  # check if it is connected
  is_connected <- igraph::is_connected(graph)
  if (!is_connected && check_connected) {
    stop("The given graph is not connected.")
  }

  if(graph_type == 'tree'){
    if(!is_tree_graph(graph)){
      stop("The given graph is not a tree.")
    }
  } else if(graph_type == 'block'){
    if(!is_block_graph(graph)){
      stop("The given graph is not a block graph.")
    }
  } else if(graph_type == 'decomposable'){
    if(!is_decomposable_graph(graph)) {
      stop("The given graph is not decomposable.")
    }
  } else if(graph_type != 'general'){
    stop("Not a valid graph_type.")
  }

  return(graph)
}

#' Check input graph and Gamma matrix
#'
#' Checks and converts the (incomplete) Gamma matrix and graph given for a
#' HR graphical model.
#'
#' @param Gamma A Gamma matrix or vector of entries corresponding to the edges
#' of `graph`
#' @param graph A graph object or `NULL` if the graph structure is specified by
#' `NA` in the Gamma matrix
#' @param graph_type Passed to [check_graph()].
#'
#' @return A list consisting of
#' \item{`Gamma`}{The Gamma matrix given as input or implied by the input}
#' \item{`graph`}{The graph given as input or implied by the input}
#' Throws an error if the input is not valid.
#'
#' @family Input checks
check_Gamma_and_graph <- function(Gamma, graph = NULL, graph_type = 'general'){
  # make graph from Gamma if necessary
  if (is.null(graph) && is.matrix(Gamma)) {
    graph <- partialMatrixToGraph(Gamma)
  } else if (is.null(graph)) {
    stop("Supply a graph or a valid Gamma matrix")
  }

  # check graph
  graph <- check_graph(graph, graph_type)

  d <- igraph::vcount(graph)
  e <- igraph::ecount(graph)

  # transform Gamma if needed
  if (is.vector(Gamma)) {
    if (length(Gamma) != e) {
      stop(paste(
        "The argument Gamma must be a symmetric d x d matrix,",
        "or a vector with as many entries as the number of edges",
        "in the graph."
      ))
    }
    G <- matrix(NA, d, d)
    edgeList <- igraph::as_edgelist(graph)
    diag(G) <- 0 # diagonal = 0
    G[edgeList] <- Gamma # upper tri
    G[edgeList[,c(2,1),drop=FALSE]] <- Gamma # lower tri
    Gamma <- G
  }

  # check that Gamma is d x d:
  if (NROW(Gamma) != d || NCOL(Gamma) != d || any(abs(Gamma - t(Gamma)) > 1e-8, na.rm = TRUE)) {
    stop(paste(
      "The argument Gamma must be a symmetric d x d matrix,",
      "or a vector with as many entries as the number of edges",
      "in the graph."
    ))
  }

  # return Gamma and graph
  return(list(
    Gamma = Gamma,
    graph = graph
  ))
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

