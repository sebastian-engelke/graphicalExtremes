
#' Completion of Gamma matrices
#'
#' Given a \code{graph} and \code{Gamma} matrix specified (at least) on the
#' edges of \code{graph}, returns the full \code{Gamma} matrix implied
#' by the conditional independencies.
#' 
#' If \code{graph} is decomposable, \code{Gamma} only needs to be specified on
#' the edges of the graph, otherwise it needs to be fully specified and a 
#' valid \{Gamma} matrix (without any graphical structure, though).
#' 
#' @param Gamma Numeric \eqn{d \times d}{d x d} variogram matrix.
#' @param graph Graph object from \code{igraph} package.
#' The \code{graph} must be a connected, undirected graph.
#' Can also be implied by \code{NA} entries in \{Gamma} if decomposable.
#' @param N Number of iterations if \code{graph} is not decomposable.
#'
#' @return Completed \eqn{d \times d}{d x d} \code{Gamma} matrix.
#'
#' @details
#' For a decomposable graph it suffices to specify the dependence parameters of the Huesler--Reiss
#' distribution within the cliques of the \code{graph}, the remaining entries are implied
#' by the conditional independence properties. For details see \insertCite{eng2019;textual}{graphicalExtremes}.
#'
#' @examples
#' ## Complete a 4-dimensional HR distribution
#'
#' my_graph <- igraph::graph_from_adjacency_matrix(rbind(
#'   c(0, 1, 0, 0),
#'   c(1, 0, 1, 1),
#'   c(0, 1, 0, 1),
#'   c(0, 1, 1, 0)
#' ),
#' mode = "undirected"
#' )
#'
#' Gamma <- rbind(
#'   c(0, .5, NA, NA),
#'   c(.5, 0, 1, 1.5),
#'   c(NA, 1, 0, .8),
#'   c(NA, 1.5, .8, 0)
#' )
#'
#' complete_Gamma(Gamma, my_graph)
#'
#' ## Alternative
#'
#' Gamma_vec <- c(.5, 1, 1.5, .8)
#' complete_Gamma(Gamma_vec, my_graph)
#' @references
#'  \insertAllCited{}
#'
#' @export
#'
complete_Gamma <- function(Gamma, graph = NULL, N = 1000, allowed_graph_type = 'general'){
  tmp <- check_Gamma_and_graph(Gamma, graph, graph_type = allowed_graph_type)
  Gamma <- tmp$Gamma
  graph <- tmp$graph
  
  if(igraph::is_chordal(graph)$chordal){
    complete_Gamma_decomposable(Gamma, graph)
  } else{
    complete_Gamma_general(Gamma, graph, N)
  }
}



complete_Gamma_general <- function(Gamma, graph, N = 1000, tol=0, check_tol=100) {

  gList <- make_graph_list(graph)
  m <- length(gList)

  for (n in 1:N) {
    t <- (n - 1) %% m + 1
    g <- gList[[t]]
    Gamma <- complete_Gamma_decomposable(Gamma, g)
    
    # Check if tolerance has been reached
    if(check_tol > 0 && n %% check_tol == 0){
      P <- Gamma2Theta(Gamma)
      A <- igraph::as_adjacency_matrix(graph)
      diag(A) <- 1
      err <- max(abs(P[A == 0]))
      if(err <= tol){
        return(Gamma)
      }
    }
  }

  return(Gamma)
}

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
  for(edge in edgeList){
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
  }
  
  return(gList)
}


complete_Gamma_decomposable <- function(Gamma, graph) {
  # computes cliques
  cliques <- igraph::max_cliques(graph)
  cliques <- order_cliques(cliques)

  # else, continue
  oldVertices <- cliques[[1]]

  # loop through remaining cliques. Skipped if only one clique in graph.
  for (p in seq_along(cliques)[-1]) {
    newVertices <- cliques[[p]]
    vC <- intersect(oldVertices, newVertices)
    vA <- setdiff(oldVertices, vC)
    vB <- setdiff(newVertices, vC)

    vACB <- c(vA, vC, vB)

    GammaACB <- Gamma[vACB, vACB]
    GammaACB <- complete_Gamma_one_step(GammaACB, length(vA), length(vC), length(vB))
    Gamma[vACB, vACB] <- GammaACB

    oldVertices <- union(oldVertices, newVertices)
  }
  return(Gamma)
}



complete_Gamma_one_step <- function(Gamma, nA, nC, nB) {
  n <- nA + nB + nC
  if (nrow(Gamma) != n || ncol(Gamma) != n) {
    stop("Make sure that nrow(Gamma) == ncol(Gamma) == nA+nB+nC")
  }

  k0 <- nA + 1 # condition on first entry in vC.
  vA <- seq_len(nA) # first nA entries
  vC <- seq_len(nC - 1) + nA # next nC-1 entries
  vB <- seq_len(nB) + nA + nC - 1 # remaining nB entries

  Sigma <- Gamma2Sigma(Gamma, k = k0)

  if (nC > 1) {
    L <- chol(Sigma[vC, vC])
    SigmaCCinv <- chol2inv(L)

    SigmaAB <- Sigma[vA, vC] %*% SigmaCCinv %*% Sigma[vC, vB]
  } else {
    SigmaAB <- 0
  }

  Sigma[vA, vB] <- SigmaAB
  Sigma[vB, vA] <- t(SigmaAB)

  Gamma <- Sigma2Gamma(Sigma, k = k0)

  return(Gamma)
}




#' Order Cliques
#'
#' Orders the cliques in a graph so that they fulfill the running intersection property.
order_cliques <- function(cliques) {
  n <- length(cliques)
  ret <- list()
  includedVertices <- numeric(0)
  for (i in 1:n) {
    for (j in seq_along(cliques)) {
      foundNextClique <- FALSE
      clique <- cliques[[j]]
      if (i == 1 || length(intersect(includedVertices, clique)) > 0) {
        # add clique to return list
        ret[[i]] <- clique
        # update list of vertices already included
        includedVertices <- union(includedVertices, clique)
        # remove clique from input list
        cliques[j] <- NULL
        foundNextClique <- TRUE
        break
      }
    }
    if (!foundNextClique) {
      stop("Graph not connected!")
    }
  }
  return(ret)
}
