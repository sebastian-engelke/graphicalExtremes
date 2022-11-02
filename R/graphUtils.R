
# Runs an igraph command (or any other expression) with
# igraph_options(add.params = FALSE)
# and resets this option afterwards
without_igraph_params <- function(expr){
  tmp <- igraph_options(add.params = FALSE)
  ret <- eval(expr, parent.frame())
  igraph_options(tmp)
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
  delete_vertex_attr(g, PERSISTENT_ID_ATTR_NAME)
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

#' Get the submatrix corresponding to a subgraph
#' 
#' Both the graph and subgraph need to have persistent IDs
#' If graph==NULL it is assumed to have pIDs 1, 2, ...
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

