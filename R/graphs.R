
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
    pids <- seq_along(ids)
  }
  g <- igraph::set_vertex_attr(g, PERSISTENT_ID_ATTR_NAME, ids, pids)
  return(g)
}
removePids <- function(g){
  delete_vertex_attr(g, PERSISTENT_ID_ATTR_NAME)
}

#' Split graph into invariant subgraphs
#' Currently cannot handle very dense graphs (e.g. complete graph on 20 nodes)
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
split_graph_at_sep <- function(g, sepIds = NULL, sepPids = NULL, returnGraphs = FALSE, includeSep = TRUE){
  if(is.null(sepIds)){
    sepIds <- getIds(g, sepPids)
  }
  sepPids <- getPids(g, sepIds)
  g2 <- igraph::delete.vertices(g, sepIds)
  compIds <- getConnectedComponents(g2)
  compsPids <- lapply(compIds, getPids, g=g2)
  if(includeSep){
    compsPids <- lapply(compsPids, c, sepPids)
  }
  partsIds <- lapply(compsPids, getIds, g=g)
  partsIds <- lapply(partsIds, sort)
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


#' Get the submatrix corresponding to a subgraph
#' 
#' The subgraph needs to have persistent IDs
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


#' Create graph list for non-decomposable completion
#' 
#' Creates a list of decomposable graphs, each consisting of two cliques,
#' such that the intersection of their edge sets
#' is identical to the edgeset of the input graph.
#' 
#' @param graph Graph object from \code{igraph} package.
#' @return List of decomposable graphs
make_graph_list <- function(graph){
  d <- igraph::vcount(graph)
  gTilde <- igraph::complementer(graph)
  edgeMat <- igraph::as_edgelist(gTilde)
  edgeList <- lapply(seq_len(NROW(edgeMat)), function(i){
    edgeMat[i,]
  })

  edgeMat0 <- igraph::as_edgelist(graph)
  edgeList0 <- lapply(seq_len(NROW(edgeMat0)), function(i){
    edgeMat0[i,]
  })

  # order edges by vertex connectivity
  conn <- sapply(edgeList, function(edge){
    igraph::vertex_connectivity(graph, edge[1], edge[2])
  })
  ind <- order(conn, decreasing = TRUE)
  edgeList <- edgeList[ind]

  gList <- list()
  partitionList <- list()
  for(edge in edgeList){
    # Check if edge is already covered by a different graph:
    alreadyCovered <- FALSE
    for(gg in gList){
      if(!igraph::are.connected(gg, edge[1], edge[2])){
        alreadyCovered <- TRUE
        break
      }
    }
    if(alreadyCovered){
      next
    }

    tmp <- igraph::max_flow(
      graph,
      edge[1],
      edge[2]
    )

    cutEdges <- edgeList0[tmp$cut]
    sep <- integer(0)
    for(ce in cutEdges){
      if(ce[1] %in% edge){
        sep <- c(sep, ce[2])
      } else {
        sep <- c(sep, ce[1])
      }
    }

    A <- union(sep, tmp$partition1)
    B <- union(sep, tmp$partition2)
    adj <- matrix(0, d, d)
    adj[A, A] <- 1
    adj[B, B] <- 1
    diag(adj) <- 0
    g <- igraph::graph_from_adjacency_matrix(adj, mode='undirected')
    gList <- c(gList, list(g))
    partitionList <- c(partitionList, list(list(
      A = A,
      B = B,
      C = intersect(A, B)
    )))
  }

  return(list(
    graphs = gList,
    partitions = partitionList
  ))
}


#' Create graph list for non-decomposable completion
#' 
#' Creates a list of decomposable graphs, each consisting of two cliques,
#' such that the intersection of their edge sets
#' is identical to the edgeset of the input graph.
#' 
#' @param graph Graph object from \code{igraph} package.
#' @return List of decomposable graphs
make_graph_list_2 <- function(graph){
  d <- igraph::vcount(graph)
  gTilde <- igraph::complementer(graph)
  edgeMat <- igraph::as_edgelist(gTilde)
  edgeList <- lapply(seq_len(NROW(edgeMat)), function(i){
    edgeMat[i,]
  })

  # get vertex connectivities
  vConns <- sapply(edgeList, function(edge){
    igraph::vertex_connectivity(graph, edge[1], edge[2])
  })
  uniqueConns <- sort(unique(vConns), decreasing = TRUE)
  
  # get minimal separators
  minSeps <- igraph::min_st_separators(graph)
  sepSizes <- sapply(minSeps, length)

  edgeIsCovered <- logical(length(edgeList))
  sepSets <- list()
  sepMats <- list()

  for(conn in uniqueConns){
    minSeps2 <- minSeps[sepSizes == conn]
    edgeList2 <- edgeList[vConns == conn & !edgeIsCovered]
    while(length(edgeList2) > 0 && length(minSeps2) > 0){
      sepMat2 <- get_sep_mat2(graph, minSeps2, edgeList2)
      ind <- which.max(rowSums(sepMat2))
      newSep <- minSeps2[[ind]]
      sepSets <- c(sepSets, list(newSep))

      edgeListFull <- edgeList[!edgeIsCovered]
      splitByNewSep <- check_split_by_sep(graph, newSep, edgeListFull)
      edgeIsCovered[!edgeIsCovered] <- splitByNewSep
      
      minSeps2 <- minSeps2[-ind]
      edgeList2 <- edgeList[vConns == conn & !edgeIsCovered]
    }
  }
  return(sepSets)
}


get_sep_mat2 <- function(graph, minSeps, edgeList){
  do.call(rbind, lapply(minSeps, check_split_by_sep, graph=graph, edgeList=edgeList))
}

check_split_by_sep <- function(graph, sep, edgeList){
  graph <- setPids(graph)
  g2 <- delete.vertices(graph, sep)
  dists <- igraph::distances(g2)
  ret <- sapply(edgeList, function(edge){
    if(any(edge %in% sep)){
      return(FALSE)
    }
    pEdge <- getPids(graph, edge)
    edge2 <- getIds(g2, pEdge)
    return(dists[edge2[1], edge2[2]] == Inf)
  })
  return(ret)
}


get_sep_mat <- function(graph, methods=c(1,2)){
  d <- igraph::vcount(graph)
  gTilde <- igraph::complementer(graph)
  edgeMat <- igraph::as_edgelist(gTilde)
  edgeList <- lapply(seq_len(NROW(edgeMat)), function(i){
    edgeMat[i,]
  })
  
  minSeps <- igraph::min_st_separators(graph)
  
  gList <- list()
  partitionList <- list()
  ms <- minSeps[[1]]
  sepMat <- matrix(NA, length(minSeps), length(edgeList))
  sepMat2 <- matrix(NA, length(minSeps), length(edgeList))
  for(i in seq_along(minSeps)){
    ms <- minSeps[[i]]
    if(1 %in% methods){
      splitGraphs <- split_graph_at_sep(g, ms, includeSep = FALSE)
      discEdges <- sapply(edgeList, is_edge_disconnected_by_partition, splitGraphs)
      sepMat[i,] <- discEdges
    }
    
    if(2 %in% methods){
      sepMat2[i,] <- check_split_by_sep(graph, ms, edgeList)
    }
  }
  
  return(list(sepMat, sepMat2))
}

is_edge_disconnected_by_partition <- function(edge, partition){
  containsEndNode <- sapply(partition, function(p){
    1 * any(edge %in% p)
  })
  return(sum(containsEndNode) == 2)
}

