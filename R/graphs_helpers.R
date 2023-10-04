
# Runs an igraph command (or any other expression) with
# igraph_options(add.params = FALSE)
# and resets this option afterwards
without_igraph_params <- function(expr){
  tmp <- igraph::igraph_options(add.params = FALSE)
  ret <- eval(expr, parent.frame())
  igraph::igraph_options(tmp)
  return(ret)
}

## Handling persistent vertex ids for graphs:
PERSISTENT_ID_ATTR_NAME <- 'pid'
getIds <- function(g, vPids = NULL){
  gPids <- igraph::vertex_attr(g, PERSISTENT_ID_ATTR_NAME)
  if(is.null(vPids)){
    return(gPids)
  }
  return(match(vPids, gPids))
}
getPids <- function(g, vIds = NULL){
  if(is.null(vIds)){
    vIds <- igraph::V(g)
  }
  return(igraph::vertex_attr(g, PERSISTENT_ID_ATTR_NAME, vIds))
}
setPids <- function(g, ids = NULL, pids = NULL, overwrite = FALSE){
  if(!overwrite && !is.null(igraph::vertex_attr(g, PERSISTENT_ID_ATTR_NAME))){
    return(g)    
  }
  if(is.null(ids)){
    ids <- igraph::V(g)
  }
  if(is.null(pids)){
    pids <- seq_along(ids)
  }
  g <- igraph::set_vertex_attr(g, PERSISTENT_ID_ATTR_NAME, ids, pids)
  return(g)
}
removePids <- function(g){
  igraph::delete_vertex_attr(g, PERSISTENT_ID_ATTR_NAME)
}


## Getting matrix entries/indices corresponding to (non-)edges in a graph
getEdgeIndices <- function(g, type = c('both', 'upper', 'lower'), withDiag = FALSE, doWhich = TRUE){
  A <- igraph::as_adjacency_matrix(g, type=type, sparse=FALSE)
  A <- (A == 1)
  if(withDiag){
    diag(A) <- TRUE
  }
  if(doWhich){
    return(which(A))
  }
  return(A)
}
getNonEdgeIndices <- function(g, type = c('both', 'upper', 'lower'), doWhich = TRUE){
  gTilde <- igraph::complementer(g)
  A <- igraph::as_adjacency_matrix(gTilde, type=type, sparse=FALSE)
  A <- (A == 1)
  if(doWhich){
    return(which(A))
  }
  return(A)
}
getEdgeEntries <- function(M, g = NULL, type = c('both', 'upper', 'lower'), withDiag = FALSE){
  if(is.null(g)){
    g <- partialMatrixToGraph(M)
  }
  A <- getEdgeIndices(g, type, withDiag)
  return(M[A])
}
getNonEdgeEntries <- function(M, g = NULL, type = c('both', 'upper', 'lower')){
  if(is.null(g)){
    g <- partialMatrixToGraph(M)
  }
  A <- getNonEdgeIndices(g, type)
  return(M[A])
}

# Creates vector of "transposed indices" in the sense
# m[ind] = t(m)[getTransposedIndices(d, ind)]
# where dim(m) == c(d, d)
getTransposedIndices <- function(d, ind = seq_len(d*d)){
  t(matrix(seq_len(d*d), d, d))[ind]
}

#' Get the submatrix corresponding to a subgraph
#' 
#' Both the graph and subgraph need to have persistent IDs
#' If graph==NULL it is assumed to have pIDs 1, 2, ...
#' @keywords internal
getSubMatrixForSubgraph <- function(fullMatrix, subgraph, graph=NULL){
  sgIds <- getIdsForSubgraph(subgraph, graph)
  return(fullMatrix[sgIds, sgIds, drop=FALSE])
}
getIdsForSubgraph <- function(subgraph, graph=NULL){
  sgPids <- getPids(subgraph)
  if(is.null(graph)){
    sgIds <- abs(sgPids) # in case negative pIDs are used
  } else{
    sgIds <- getIds(graph, sgPids)
  }
  return(sgIds)
}

## Check if a graph is of a certain type
is_block_graph <- function(graph, check_connected=TRUE){
  if(check_connected && !igraph::is_connected(graph)){
    return(FALSE)
  }
  if(!is_decomposable_graph(graph)){
    return(FALSE)
  }
  # Check that separators are size 1 or 0:
  cliques <- igraph::max_cliques(graph)
  for(i in seq_along(cliques)){
    for(j in seq_len(i-1)){
      if(length(intersect(cliques[[i]], cliques[[j]])) > 1){
        return(FALSE)
      }
    }
  }
  return(TRUE)
}
is_tree_graph <- function(graph){
  igraph::is_tree(graph, 'all', details = FALSE)
}
is_decomposable_graph <- function(graph){
  igraph::is.chordal(graph)$chordal
}
is_complete_graph <- function(graph){
  d <- igraph::vcount(graph)
  return(igraph::ecount(graph) == d*(d-1)/2)
}

#' Graph equality
#' 
#' Produce true if two graphs have same vertices and edges (ordered)
#' 
#' @param g1 `igraph::graph`
#' @param g2 `igraph::graph`
#' @return `logical` indicating if the graphs are equal
#' @export
graphs_equal <- function(g1, g2) {

  # Return early if graph sizes are different
  if(igraph::vcount(g1) != igraph::vcount(g2)){
    return(FALSE)
  }
  # Compare adjacency matrices
  A1 <- igraph::as_adjacency_matrix(g1, sparse = FALSE)
  A2 <- igraph::as_adjacency_matrix(g2, sparse = FALSE)
  return(all(A1 == A2))
}
