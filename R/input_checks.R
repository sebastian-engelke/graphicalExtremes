#' Check input graph
#'
#' Checks that the input graph is a valid graph for an extremal graphical model.
#' If necessary, converts the graph to an undirected graph.
#'
#' @param graph An [igraph::graph] object
#' @param graph_type `"general"`, `"decomposable"`, `"block"`, `"tree"`
#' The required type of graph
#' @param check_connected Whether to check if the graph is connected
#' @param nVertcies The number of vertices required in the graph
#'
#' @return The given `graph`, if necessary converted to undirected.
#' If the graph is not valid an error is thrown.
#'
#' @family Input checks
check_graph <- function(graph, graph_type='general', check_connected=TRUE, nVertices=NULL){

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
    warning("The given graph is directed. Converted to undirected.")
    graph <- igraph::as.undirected(graph)
  }

  if(!is.null(igraph::vertex_attr(graph)[['name']])){
    warning("The vertex labels were removed.")
    igraph::vertex_attr(graph)[['name']] <- NULL
  }

  # check if it is connected
  is_connected <- igraph::is_connected(graph)
  if (!is_connected && check_connected) {
    stop("The given graph is not connected.")
  }

  is_tree <- is_connected && e == d-1
  is_block <- is_block_graph(graph)
  is_decomposable <- igraph::is_chordal(graph)$chordal
  if(graph_type == 'tree'){
    if(!is_tree){
      stop("The given graph is not a tree.")
    }
  } else if(graph_type == 'block'){
    if(!is_block){
      stop("The given graph is not a block graph.")
    }
  } else if(graph_type == 'decomposable'){
    if(!is_decomposable) {
      stop("The given graph is not decomposable.")
    }
  } else if(graph_type != 'general'){
    stop("Not a valid graph_type.")
  }

  return(graph)
}

#' Check if a graph is a block graph
#'
#' @param graph An [igraph::graph] object
#'
#' @return A `boolean` indicating if the graph is a glock graph
#'
#' @family Input checks
is_block_graph <- function(graph){
  if(!igraph::is_connected(graph) || !igraph::is_chordal(graph)$chordal){
    return(FALSE)
  }
  cliques <- igraph::max_cliques(graph)

  # Check that separators are size 1 or 0:
  for(i in seq_along(cliques)){
    for(j in seq_len(i-1)){
      if(length(intersect(cliques[[i]], cliques[[j]])) > 1){
        return(FALSE)
      }
    }
  }
  return(TRUE)
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
#' @return A list containing two named entries, `Gamma` and `graph` containing
#' the input matrix and graph. Throws an error if the input is not valid.
#'
#' @family Input checks
check_Gamma_and_graph <- function(Gamma, graph = NULL, graph_type = 'general'){

  # make graph from Gamma if necessary
  if (is.null(graph) && is.matrix(Gamma)) {
    graph <- igraph::graph_from_adjacency_matrix(
      1 * !is.na(Gamma),
      mode = "undirected"
    )
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
    G <- matrix(0, d, d)
    G[igraph::as_edgelist(graph)] <- Gamma
    Gamma <- G + t(G)
  }

  # check that Gamma is d x d:
  if (NROW(Gamma) != d || NCOL(Gamma) != d || any(abs(Gamma - t(Gamma)) > 1e-8, na.rm = T)) {
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


