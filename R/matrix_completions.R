
#' Completion of Gamma matrices
#'
#' Given a `graph` and `Gamma` matrix specified (at least) on the
#' edges of `graph`, returns the full `Gamma` matrix implied
#' by the conditional independencies.
#'
#' If `graph` is decomposable, `Gamma` only needs to be specified on
#' the edges of the graph and the graph structure can be implied by setting
#' the remaining entries to `NA`.
#'
#' If `graph` is not decomposable, the algorithm requires a fully specified
#' variogram matrix `Gamma` and the graph structure needs to be explicitly
#' provided in `graph`.
#'
#' @param Gamma Numeric \eqn{d \times d}{d x d} variogram matrix.
#' @param graph Graph object from \code{igraph} package.
#' The \code{graph} must be a connected, undirected graph.
#' Can also be implied by \code{NA} entries in \code{Gamma} if decomposable.
#' @param allowed_graph_type Is passed as `graph_type` to [check_graph()].
#' Can be used to throw an error if `graph` is not of the specified type,
#' but does not have any influence on the completion algorithm.
#' @param ... Further arguments passed to [complete_gamma_general()] if `graph`
#' is not decomposable
#'
#' @return Completed \eqn{d \times d}{d x d} \code{Gamma} matrix.
#'
#' @details
#' For a decomposable graph it suffices to specify the dependence parameters of the Huesler--Reiss
#' distribution within the cliques of the \code{graph}, the remaining entries are implied
#' by the conditional independence properties. For details see \insertCite{eng2019;textual}{graphicalExtremes}.
#'
#' @examples
#' ## Block graph:
#' Gamma <- rbind(
#'   c(0, .5, NA, NA),
#'   c(.5, 0, 1, 1.5),
#'   c(NA, 1, 0, .8),
#'   c(NA, 1.5, .8, 0)
#' )
#'
#' complete_Gamma(Gamma)
#'
#' ## Alternative representation of the same completion problem:
#' my_graph <- igraph::graph_from_adjacency_matrix(rbind(
#'   c(0, 1, 0, 0),
#'   c(1, 0, 1, 1),
#'   c(0, 1, 0, 1),
#'   c(0, 1, 1, 0)
#' ), mode = "undirected")
#' Gamma_vec <- c(.5, 1, 1.5, .8)
#' complete_Gamma(Gamma_vec, my_graph)
#'
#' ## Decomposable graph:
#' G <- rbind(
#' c(0, 5, 7, 6, NA),
#' c(5, 0, 14, 15, NA),
#' c(7, 14, 0, 5, 5),
#' c(6, 15, 5, 0, 6),
#' c(NA, NA, 5, 6, 0)
#' )
#'
#' complete_Gamma(G)
#'
#' ## Non-decomposable graph:
#' G <- rbind(
#' c(0, 5, 7, 6, 6),
#' c(5, 0, 14, 15, 13),
#' c(7, 14, 0, 5, 5),
#' c(6, 15, 5, 0, 6),
#' c(6, 13, 5, 6, 0)
#' )
#' g <- igraph::make_ring(5)
#'
#' complete_Gamma(G, g)
#'
#'
#' @references
#'  \insertAllCited{}
#'
#' @seealso [Gamma2Theta()]
#' @family Matrix completions
#' @export
#'
complete_Gamma <- function(Gamma, graph = NULL, allowed_graph_type = 'general', ...){
  tmp <- check_Gamma_and_graph(Gamma, graph, graph_type = allowed_graph_type)
  Gamma <- tmp$Gamma
  graph <- tmp$graph

  if(igraph::is_chordal(graph)$chordal){
    complete_Gamma_decomposable(Gamma, graph)
  } else{
    complete_Gamma_general(Gamma, graph, ...)
  }
}


#' Completion of non-decomposable Gamma matrices (SLOW)
#'
#' Given a \code{graph} and variogram matrix `Gamma`, returns the full \code{Gamma}
#' matrix implied by the conditional independencies.
#' This function uses a convergent iterative algorithm.
#' SLOW VERSION. Allows returning details, specifying graph list.
#'
#' @param Gamma A complete variogram matrix (without any graphical structure)
#' @param graph An [igraph::graph] object
#' @param N The maximal number of iterations of the algorithm
#' @param tol The tolerance to use when checking for zero entries in `Theta`
#' @param check_tol After how many iterations to check the tolerance in `Theta`
#'
#' @return A matrix that agrees with `Gamma` on the entries corresponding to
#' edges in `graph` and the diagonals.
#' The corresponding \eqn{\Theta} matrix produced by [Gamma2Theta] has values
#' close to zero in the remaining entries (how close depends on the input
#' and the number of iterations).
#'
#' @family Matrix completions
complete_Gamma_general_0 <- function(Gamma, graph, N = 1000, tol=0, check_tol=100, saveDetails=FALSE, gList=NULL) {

  if(is.null(gList)){
    gList <- make_graph_list(graph)$graphs
  }
  m <- length(gList)
  
  if(saveDetails){
    GammaList <- list(Gamma)
  }

  for (n in 1:N) {
    t <- (n - 1) %% m + 1
    g <- gList[[t]]
    Gamma <- complete_Gamma_decomposable(Gamma, g)
    
    if(saveDetails){
      GammaList <- c(GammaList, list(Gamma))
    }

    # Check if tolerance has been reached
    if(check_tol > 0 && n %% check_tol == 0){
      P <- Gamma2Theta(Gamma)
      A <- igraph::as_adjacency_matrix(graph, sparse=FALSE)
      diag(A) <- 1
      err <- max(abs(P[A == 0]))
      if(err <= tol){
        break
      }
    }
  }

  if(saveDetails){
    return(list(
      GammaList = GammaList,
      graphList = gList,
      N = N,
      tol = tol,
      check_tol = check_tol
    ))
  }

  return(Gamma)
}



#' Completion of non-decomposable Gamma matrices (FASTER)
#'
#' Given a \code{graph} and variogram matrix `Gamma`, returns the full \code{Gamma}
#' matrix implied by the conditional independencies.
#' This function uses a convergent iterative algorithm.
#'
#' @param Gamma A complete variogram matrix (without any graphical structure)
#' @param graph An [igraph::graph] object
#' @param N The maximal number of iterations of the algorithm
#' @param tol The tolerance to use when checking for zero entries in `Theta`
#' @param check_tol After how many iterations to check the tolerance in `Theta`
#'
#' @return A matrix that agrees with `Gamma` on the entries corresponding to
#' edges in `graph` and the diagonals.
#' The corresponding \eqn{\Theta} matrix produced by [Gamma2Theta] has values
#' close to zero in the remaining entries (how close depends on the input
#' and the number of iterations).
#'
#' @family Matrix completions
#' @export
complete_Gamma_general <- function(Gamma, graph, N = 1000, tol=0, check_tol=100) {

  tmp <- make_graph_list(graph)
  partitionList <- tmp$partitions
  gList <- tmp$graphs
  
  indList <- lapply(partitionList, function(AB) list(
    vC = intersect(AB$A, AB$B),
    vA = setdiff(AB$A, AB$B),
    vB = setdiff(AB$B, AB$A)
  ))
  m <- length(indList)
  
  if(length(gList) == 0){
    N <- 0
  }

  for (n in seq_len(N)) {
    t <- (n - 1) %% m + 1
    vABC <- indList[[t]]
    vA <- vABC$vA
    vB <- vABC$vB
    vC <- vABC$vC
    
    if(length(vC) == 1){
      k0 <- vC[1]
      GammaAB <- outer(Gamma[vA, k0], Gamma[k0, vB], '+')
      Gamma[vA, vB] <- GammaAB
      Gamma[vB, vA] <- t(GammaAB)
    } else{
      k0 <- vC[1]
      vC_Sigma <- vC[-1]
      Sigma <- Gamma2Sigma(Gamma, k = k0, full = TRUE)
      R <- chol(Sigma[vC_Sigma, vC_Sigma, drop=FALSE])
      SigmaCCinv <- chol2inv(R)
      SigmaAB <- Sigma[vA, vC_Sigma, drop=FALSE] %*% SigmaCCinv %*% Sigma[vC_Sigma, vB, drop=FALSE]
      Sigma[vA, vB] <- SigmaAB
      Sigma[vB, vA] <- t(SigmaAB)
      Gamma <- Sigma2Gamma(Sigma)
    }
    
    # Check if tolerance has been reached
    if(check_tol > 0 && n %% check_tol == 0){
      P <- Gamma2Theta(Gamma)
      A <- igraph::as_adjacency_matrix(graph, sparse=FALSE)
      diag(A) <- 1
      err <- max(abs(P[A == 0]))
      if(err <= tol){
        break
      }
    }
  }

  return(Gamma)
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
      B = B
    )))
  }

  return(list(
    graphs = gList,
    partitions = partitionList
  ))
}

#' Completion of decomposable Gamma matrices
#'
#' Given a decomposable `graph` and incomplete variogram matrix `Gamma`,
#' returns the full `Gamma` matrix implied by the conditional independencies.
#'
#' @param Gamma A variogram matrix that is specified on the edges of `graph`
#' and the diagonals. All other entries are ignored.
#' @param graph A decomposable [igraph::graph] object
#'
#' @return A complete variogram matrix that agrees with `Gamma` on the entries
#' corresponding to edges in `graph` and the diagonals.
#' The corresponding \eqn{\Theta} matrix pdocued by [Gamma2Theta] has zeros
#' in the remaining entries.
#'
#' @family Matrix completions
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



#' Completion of two-clique decomposable Gamma matrices
#'
#' Given a decomposable `graph` consisting of two cliques and incomplete
#' variogram matrix `Gamma`,
#' returns the full `Gamma` matrix implied by the conditional independencies.
#'
#' @param Gamma A variogram matrix that is specified on the edges of `graph`
#' and the diagonals. All other entries are ignored.
#' @param graph A decomposable [igraph::graph] object with two cliques and
#' non-empty separator set
#'
#' @return A complete variogram matrix that agrees with `Gamma` on the entries
#' corresponding to edges in `graph` and the diagonals.
#' The corresponding \eqn{\Theta} matrix pdocued by [Gamma2Theta] has zeros
#' in the remaining entries.
#'
#' @family Matrix completions
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
    R <- chol(Sigma[vC, vC, drop=FALSE])
    SigmaCCinv <- chol2inv(R)
    SigmaAB <- Sigma[vA, vC] %*% SigmaCCinv %*% Sigma[vC, vB]
    
    ## These don't seem to improve performance:
    # # X <- solve(Sigma[vC, vC, drop=FALSE], Sigma[vC, vB, drop=FALSE])
    # Y <- backsolve(R, Sigma[vC, vB, drop=FALSE])
    # X <- forwardsolve(t(R), Y)
    # SigmaAB <- Sigma[vA, vC, drop=FALSE] %*% X
  } else {
    SigmaAB <- 0
  }

  Sigma[vA, vB] <- SigmaAB
  Sigma[vB, vA] <- t(SigmaAB)

  Gamma <- Sigma2Gamma(Sigma, k = k0)

  return(Gamma)
}





try_complete_Gamma <- function(graph, Gamma, key, val){
  ## igraph numeric_matrix character double -> numeric_matrix | NA
  ## tries to call `complete_Gamma`, if it fails returns NA

  # Try to complete Gamma
  Gamma_curr <- tryCatch(
    {complete_Gamma(graph = graph, Gamma = Gamma)},
    error = function(e){
      if (e$message == "The given graph is not connected.") {
        message(glue::glue(
          "The estimated graph for {key} = {val}",
          " is not connected,",
          " so it is not possible to complete `Gamma`.\n"
        ))
        NA
      }
      else {
        stop(e)
      }
    }
  )

  # Check if completed Gamma matches with given graph
  if (all(!is.na(Gamma_curr))) {
    completed_graph <- Gamma2graph(Gamma_curr, to_plot = FALSE)
    if (!(graphs_equal(completed_graph, graph))) {
      message(glue::glue(
        "The completed Gamma for {key} = {val}",
        " does not match the estimated graph.\n"
      ))
    }
  }

  # Return completed Gamma
  return(Gamma_curr)
}
