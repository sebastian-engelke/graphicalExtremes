
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
#' @param graph Graph object from `igraph` package.
#' The `graph` must be a connected, undirected graph.
#' Can also be implied by `NA` entries in `Gamma` if decomposable.
#' @param allowed_graph_type Is passed as `graph_type` to [check_graph()].
#' Can be used to throw an error if `graph` is not of the specified type,
#' but does not have any influence on the completion algorithm.
#' @param ... Further arguments passed to [complete_gamma_general()] if `graph`
#' is not decomposable
#'
#' @return Completed \eqn{d \times d}{d x d} `Gamma` matrix.
#'
#' @details
#' For a decomposable graph it suffices to specify the dependence parameters of the Huesler--Reiss
#' distribution within the cliques of the `graph`, the remaining entries are implied
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
complete_Gamma <- function(
  Gamma,
  graph = NULL,
  allowed_graph_type = 'general',
  ...
){
  tmp <- check_Gamma_and_graph(Gamma, graph, graph_type = allowed_graph_type)
  Gamma <- tmp$Gamma
  graph <- tmp$graph

  # Return completion if graph is decomposable (=chordal)
  if(igraph::is_chordal(graph)$chordal){
    return(complete_Gamma_decomposable(Gamma, graph))
  }

  # Compute initial non-graphical completion if necessary:
  if(any(is.na(Gamma))){
    A <- 1*!is.na(Gamma)
    tmp <- edmcr::npf(Gamma, A, d = NROW(Gamma)-1)
    Gamma <- tmp$D
    if(!is_sym_cnd(Gamma)){
      stop('Did not find an initial non-graphical completion (using edmcr::npf)!')
    }
  }

  # Compute non-decomposable completion
  return(complete_Gamma_general_mc(Gamma, graph, ...))
}



#' Completion of non-decomposable Gamma matrices (demo-version)
#'
#' Given a `graph` and variogram matrix `Gamma`, returns the full `Gamma`
#' matrix implied by the conditional independencies.
#' This function uses a convergent iterative algorithm.
#' DEMO VERSION: Returns a lot of details and allows specifying the graph list
#' that is used. Is way slower than `complete_Gamma_general`.
#'
#' @param Gamma A complete variogram matrix (without any graphical structure)
#' @param graph An [igraph::graph] object
#' @param N The maximal number of iterations of the algorithm
#' @param tol The tolerance to use when checking for zero entries in `Theta`
#'
#' @return A nested list, containing the following details: TODO!!!
#' 
#' The corresponding \eqn{\Theta} matrix produced by [Gamma2Theta] has values
#' close to zero in the remaining entries (how close depends on the input
#' and the number of iterations).
#'
#' @family Matrix completions
complete_Gamma_general_demo <- function(Gamma, graph = NULL, N = 1000, tol=0, gList=NULL) {
  # Compute gList if not provided:
  if(is.null(gList)){
    sepDetails <- make_sep_list(graph, details=TRUE)
    gList <- lapply(sepDetails, function(dets) dets$graph)
  }
  m <- length(gList)

  # Initialize ret-list with initial Gamma, Theta, etc.:
  Theta <- Gamma2Theta(Gamma)
  iterations <- list()
  ret <- list(
    graph = graph,
    Gamma0 = Gamma,
    Theta0 = Theta,
    err0 = max(abs(getNonEdgeEntries(Theta, graph))),
    gList = gList,
    tol = tol,
    N = N
  )

  # Iterate over gList:
  n <- 0
  while(n < N) {
    n <- n + 1
    t <- (n - 1) %% m + 1
    g <- gList[[t]]
    Gamma <- complete_Gamma_decomposable(Gamma, g)

    GammaList <- c(GammaList, list(Gamma))

    # Compute Theta
    Theta <- Gamma2Theta(Gamma)
    err <- max(abs(getNonEdgeEntries(Theta, graph)))

    # Store results
    iterations <- c(iterations, list(list(
      n = n,
      t = t,
      g = g,
      Gamma = Gamma,
      Theta = Theta,
      err = err
    )))

    # Check if tolerance has been reached
    if(err <= tol){
      break
    }
  }
  
  ret$iterations <- iterations

  return(ret)
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
complete_Gamma_decomposable <- function(Gamma, graph = NULL) {
  # Compute graph if not specified
  if(is.null(graph)){
    graph <- partialMatrixToGraph(Gamma)
  }

  # compute cliques and order by the running intersection property
  cliques <- igraph::max_cliques(graph)
  cliques <- order_cliques(cliques)

  # Start with first clique
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
  Gamma <- ensure_symmetry(Gamma) # todo: set tolerance = Inf?
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
  # Check arguments
  n <- nA + nB + nC
  if (nrow(Gamma) != n || ncol(Gamma) != n) {
    stop("Make sure that nrow(Gamma) == ncol(Gamma) == nA+nB+nC")
  }

  # Prepare index sets
  vA <- seq_len(nA) # first nA entries
  vC <- seq_len(nC) + nA # next nC entries
  vB <- seq_len(nB) + nA + nC # remaining nB entries
  k0 <- vC[1] # condition on first entry in vC.
  vC_Sigma <- vC[-1]
  
  if(nC == 1){
    ## Separator of size 1 -> use additive property of block matrix completion
    GammaAB <- outer(Gamma[vA, k0], Gamma[k0, vB], `+`)
    Gamma[vA, vB] <- GammaAB
    Gamma[vB, vA] <- t(GammaAB)
    return(Gamma)
  }

  ## Larger separator -> convert to Sigma, complete, convert back
  Sigma <- Gamma2Sigma(Gamma, k = k0, full = TRUE)

  # Invert separator submatrix
  R <- chol(Sigma[vC_Sigma, vC_Sigma, drop=FALSE])
  SigmaCCinv <- chol2inv(R)
  # Compute completion
  SigmaAB <- Sigma[vA, vC_Sigma] %*% SigmaCCinv %*% Sigma[vC_Sigma, vB]

  # Fill in completed values
  Sigma[vA, vB] <- SigmaAB
  Sigma[vB, vA] <- t(SigmaAB)

  # Convert back
  Gamma <- Sigma2Gamma(Sigma)
  return(Gamma)
}



