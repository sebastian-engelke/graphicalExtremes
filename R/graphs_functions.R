
## This file contains functions related (solely) to graph manipulation

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


#' Get Cliques and Separators of a graph
#'
#' Finds all cliques, separators, and (recursively) separators of separators
#' in a graph.
#'
#' @param graph An [igraph::graph] object
#' @return A list of vertex sets that represent the cliques and (recursive)
#' separators of `graph`
#'
#' @keywords internal
get_cliques_and_separators_OLD <- function(graph){
  # start with maximal cliques as new cliques:
  newCliques <- igraph::max_cliques(graph)
  cliques <- list()

  while(length(newCliques)>0){
    # compute all separators between the new cliques:
    separators <- get_separators(newCliques)
    # add new cliques to cliques-list:
    # (separators are added after the next iteration)
    cliques <- c(newCliques, cliques)

    # use separators as new cliques:
    newCliques <- setdiff(separators, cliques)
  }

  return(cliques)
}


#' Get Cliques and Separators of a graph
#'
#' Finds all cliques, separators, and (recursively) separators of separators
#' in a graph.
#'
#' @param graph An [igraph::graph] object
#' @return A list of vertex sets that represent the cliques and (recursive)
#' separators of `graph`, ordered such that separators come before cliques they separate.
#'
#' @keywords internal
get_cliques_and_separators <- function(graph, sortIntoLayers = FALSE, includeSingletons = FALSE){
  # start with maximal cliques of graph
  newCliques <- lapply(igraph::max_cliques(graph), as.numeric)
  allCliques <- newCliques

  while(length(newCliques) > 0){
    # compute all separators
    separators <- get_separators(allCliques, includeSingletons)
    # check which are actually new cliques
    newCliques <- setdiff(separators, allCliques)
    # add to list of cliques
    allCliques <- c(newCliques, allCliques)
  }
  if(!sortIntoLayers){
    return(allCliques)
  }
  layers <- sort_cliques_and_separators(allCliques)
  return(layers)
}

#' Sort a set of cliques and separators into layers, such that the separator
#' of two cliques/separators is contained in a previous layer
sort_cliques_and_separators <- function(cliques){
  # Compute matrix indicating whether clique[[i]] is a subset of clique[[j]]
  subsetMat <- matrix(FALSE, length(cliques), length(cliques))
  for(i in seq_along(cliques)){
    for(j in seq_along(cliques)){
      if(i != j){
        subsetMat[i,j] <- all(cliques[[i]] %in% cliques[[j]])
      }
    }
  }
  # Make layers
  layers <- list()
  while(nrow(subsetMat) > 0){
    isSmallest <- apply(subsetMat, 2, function(v) !any(v))
    if(!any(isSmallest)){
      stop('hae?')
    }
    layers <- c(layers, list(cliques[isSmallest]))
    subsetMat <- subsetMat[!isSmallest, !isSmallest, drop=FALSE]
    cliques <- cliques[!isSmallest]
  }
  return(layers)
}

# Get all separators (non-recursive) between a set of cliques
get_separators <- function(cliques, includeSingletons = FALSE){
  separators <- list()
  for(i in seq_along(cliques)){
    for(j in seq_len(i-1)){
      sep <- sort(intersect(cliques[[i]], cliques[[j]]))
      if(length(sep)>1 || (includeSingletons && length(sep) == 1)){
        separators <- c(separators, list(sep))
      }
    }
  }
  separators <- unique(separators)
  return(separators)
}


#' Split graph into invariant subgraphs
split_graph <- function(g){
  g <- setPids(g)
  # compute cliques and sort by size:
  cliques <- igraph::max_cliques(g)
  ord <- order(sapply(cliques, length))
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
      if(!igraph::is.separator(g, cli)){
        newGraphs <- c(newGraphs, list(g))
        next # cli not a separator in g
      }
      splitGraphs <- split_graph_at_sep(g, cli, returnGraphs = TRUE)
      newGraphs <- c(newGraphs, splitGraphs)
    }
    graphs <- newGraphs
  }
  graphs <- split_off_cliques(graphs, pCliques)
  return(graphs)
}

split_off_cliques <- function(graphs, pCliques){
  for(pCli in pCliques){
    newGraphs <- list()
    for(g in graphs){
      if(!all(pCli %in% getPids(g))){
        newGraphs <- c(newGraphs, list(g))
        next # cli not fully contained in g
      }
      cli <- getIds(g, pCli)
      cliDegrees <- igraph::degree(g, cli)
      ind <- (cliDegrees > length(cli) - 1) # vertices in cli, connected to the rest of g
      if(all(ind)){
        newGraphs <- c(newGraphs, list(g))
        next # all vertices of cli are connected to the rest of g
      }
      newSep <- cli[ind]
      splitGraphs <- split_graph_at_sep(g, newSep, returnGraphs = TRUE)
      newGraphs <- c(newGraphs, splitGraphs)
    }
    graphs <- newGraphs
  }

  # Order graphs large -> small
  graphs <- graphs[order(sapply(graphs, igraph::vcount), decreasing = TRUE)]

  # Remove cliques that are complete subsets (or equal) of another graph:
  newGraphs <- list()
  for(g in graphs){
    contains_g <- sapply(newGraphs, function(g0) {
      all(getPids(g) %in% getPids(g0))
    })
    if(!any(contains_g)){
      newGraphs <- c(newGraphs, list(g))
    }
  }
  return(newGraphs)
}

#' Split a graph at a single separator
split_graph_at_sep <- function(g, sepIds = NULL, sepPids = NULL, returnGraphs = FALSE, includeSep = TRUE){
  g <- setPids(g)
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


#' Create graph list for non-decomposable completion
#' 
#' Creates a list of decomposable graphs, each consisting of two cliques,
#' such that the intersection of their edge sets
#' is identical to the edgeset of the input graph.
#' 
#' @param graph Graph object from `igraph` package.
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


#' Create a list of separators
#' 
#' Creates a list of separator set, such that every pair of non-adjacent
#' vertices in `graph` is completely disconnected by the removal of
#' (at least) one of the separator sets from the graph.
#' 
#' @param graph A graph
#' @return A list of numeric vectors 
make_sep_list <- function(graph, details=TRUE){
  graph <- setPids(graph)
  # Get matrix of non-edges (edgeTildes):
  gTilde <- igraph::complementer(graph)
  edgeTildeMat <- igraph::as_edgelist(gTilde)

  # order edgeTildes by vertex connectivity (in the original graph):
  conn <- sapply(seq_len(NROW(edgeTildeMat)), function(i){
    igraph::vertex_connectivity(graph, edgeTildeMat[i,1], edgeTildeMat[i,2])
  })
  ind <- order(conn, decreasing = TRUE)
  edgeTildeMat <- edgeTildeMat[ind,]

  # All edgeTildes need to be separated (=covered) by some sep:
  edgeTildeIsCovered <- logical(NROW(edgeTildeMat))
  sepList <- list()

  # Loop over edgeTildes, find separators:
  while(any(!edgeTildeIsCovered)){
    # Find first edgeTilde that is not covered yet:
    i <- match(FALSE, edgeTildeIsCovered)

    # Find a separator set for this edgeTilde:
    vSep <- findVsep(graph, edgeTildeMat[i,1], edgeTildeMat[i,2])
    edgeTildeIsCovered[i] <- TRUE
    sepList <- c(sepList, list(vSep))

    # Find all other edgeTildes that are split by this sep:
    splitEdgeTildes <- check_split_by_sep(
      graph,
      vSep,
      edgeTildeMat[!edgeTildeIsCovered,,drop=FALSE]
    )
    edgeTildeIsCovered[!edgeTildeIsCovered] <- splitEdgeTildes
  }

  if(!details){
    return(sepList)
  }
  # Add details for each sep:
  detailedSepList <- lapply(sepList, makeSepDetails, graph=graph)
  return(detailedSepList)
}

makeSepDetails <- function(graph, sep){
  parts <- split_graph_at_sep(graph, sep, includeSep = FALSE)
  partsWithSep <- lapply(parts, function(part){
    sort(c(part, sep))
  })
  partPairs <- utils::combn(parts, 2, simplify = FALSE)
  k <- sep[1]
  sepWithoutK <- sep[-1]
  d <- igraph::vcount(graph)
  A <- matrix(0, d, d)
  for(part in partsWithSep){
    A[part, part] <- 1
  }
  g2 <- igraph::graph_from_adjacency_matrix(A, 'undirected', diag = FALSE)
  return(list(
    parts = parts,
    partsWithSep = partsWithSep,
    partPairs = partPairs,
    k = k,
    sep = sep,
    sepWithoutK = sepWithoutK,
    graph = g2
  ))
}

#' Find a separator set for two vertices
#' 
#' Finds a reasonably small set of vertices that separate `v0` and `v1` in `graph`.
findVsep <- function(graph, v0, v1){
  # find max flow between vertices (includes edge-sep):
  tmp <- igraph::max_flow(
    graph,
    v0,
    v1
  )

  # convert edge-sep to vertex-sep:
  edgeMat <- igraph::as_edgelist(graph)
  eSep <- edgeMat[tmp$cut,,drop=FALSE]
  useSecondVertex <- (eSep[,1] %in% c(v0, v1))
  vSep <- eSep[,1]
  vSep[useSecondVertex] <- eSep[useSecondVertex,2]
  
  vSep <- unique(vSep)

  return(vSep)
}

#' Identify pairs of vertices that are split by a separator
#' 
#' @param graph A graph
#' @param sep A set of vertex ids that are used to split the graph
#' @param edgeMat A two-column matrix, containing the vertex-paris to be checked
#' 
#' @return A logical vector, indicating for each edge whether it is split by `sep`
#' 
check_split_by_sep <- function(graph, sep, edgeMat){
  # Return early if no vertex pair (edge) is to be checked:
  if(nrow(edgeMat) == 0){
    return(logical(0))
  }

  # Compute distances after removing `sep`:
  g2 <- igraph::delete.vertices(graph, sep)
  dists <- igraph::distances(g2)

  # Convert vertex ids in original graph to ids in split graph:
  graph <- setPids(graph)
  pEdgeMat <- matrix(getPids(graph, edgeMat), ncol = 2)
  edgeMat2 <- matrix(getIds(g2, pEdgeMat), ncol = 2)

  # Read distances between each vertex pair in the split graph:
  edgeDists <- dists[edgeMat2]
  edgeVertexInSep <- matrix(edgeMat %in% sep, ncol = 2)

  # Vertexes are split if their distance is Inf in the split graph:
  isSplit <- (edgeDists == Inf) & !is.na(edgeDists)
  return(isSplit)
}

