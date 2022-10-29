
## This file contains functions related (solely) to graph manipulation

PERSISTENT_ID_ATTR_NAME <- 'pid'

#' Order Cliques
#'
#' Orders the cliques in a connected decomposable graph so that they fulfill the running intersection property.
order_cliques <- function(cliques) {
  n <- length(cliques)
  ret <- list()
  for (i in 1:n) {
    foundNextClique <- FALSE
    for (j in seq_along(cliques)) {
      candidate <- cliques[[j]]
      rest <- Reduce(union, cliques[-j], c())
      separator <- intersect(candidate, rest)
      sepInClique <- sapply(cliques[-j], function(C) all(separator %in% C))
      if (length(sepInClique) == 0 || any(sepInClique)) {
        # add clique to return list
        # ret[[i]] <- candidate
        ret <- c(list(candidate), ret)
        # remove clique from input list
        cliques[j] <- NULL
        foundNextClique <- TRUE
        break
      }
    }
    if (!foundNextClique) {
      stop("Graph not decomposable or not connected!")
    }
  }
  return(ret)
}


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
    pids <- (-1) * seq_along(ids) # make negative to distinguish easily from ids
  }
  g <- igraph::set_vertex_attr(g, PERSISTENT_ID_ATTR_NAME, ids, pids)
  return(g)
}
removePids <- function(g){
  delete_vertex_attr(g, PERSISTENT_ID_ATTR_NAME)
}

#' Split graph into invariant subgraphs
split_graph <- function(g){
  g <- setPids(g)
  # compute cliques and sort by size:
  cliques <- igraph::cliques(g)
  ord <- order(sapply(cliques(g), length))
  cliques <- cliques[ord]
  pCliques <- lapply(cliques, getPids, g=g)

  graphs <- list(g)
  for(pCli in pCliques){
    newGraphs <- list()
    for(g in graphs){
      if(!all(pCli %in% getPids(g))){
        newGraphs <- c(newGraphs, list(g))
        next # cli not fully contained in g
      }
      cli <- getIds(g, pCli)
      if(!is.separator(g, cli)){
        newGraphs <- c(newGraphs, list(g))
        next # cli not a separator in g
      }
      splitGraphs <- split_graph_at_sep(g, cli, returnGraphs = TRUE)
      newGraphs <- c(newGraphs, splitGraphs)
    }
    graphs <- newGraphs
  }
  return(graphs)
}

#' Split a graph at a single separator
split_graph_at_sep <- function(g, sepIds = NULL, sepPids = NULL, returnGraphs = FALSE){
  if(is.null(sepIds)){
    sepIds <- getIds(g, sepPids)
  }
  sepPids <- getPids(g, sepIds)
  g2 <- igraph::delete.vertices(g, sepIds)
  compIds <- getConnectedComponents(g2)
  compsPids <- lapply(compIds, getPids, g=g2)
  partsPids <- lapply(compsPids, function(compPids){
    sort(c(sepPids, compPids))
  })
  partsIds <- lapply(partsPids, getIds, g=g)
  if(returnGraphs){
    return(lapply(partsIds, igraph::induced_subgraph, graph=g))
  }
  return(partsIds)
}

getConnectedComponents <- function(g){
  tmp <- igraph::components(g)
  components <- lapply(seq_len(tmp$no), function(i){
    which(tmp$membership == i)
  })
  return(components)
}



